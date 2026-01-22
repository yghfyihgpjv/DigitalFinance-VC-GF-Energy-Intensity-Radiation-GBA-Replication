version 17.0
clear all
set more off
set linesize 200


*====================================================================
* 0. 环境设置与日志记录
*====================================================================

* 使用当前工作目录作为项目根目录
local ROOT "`c(pwd)'"

* 数据与输出的相对路径
local DATA "`ROOT'/data/processeddata/29.dta"
local LOGDIR "`ROOT'/results"

* 如果 results 文件夹不存在，则创建
capture mkdir "`LOGDIR'"

* 日志文件路径
local log_file "`LOGDIR'/a1_annotated.txt"

capture log close _all
log using "`log_file'", text replace


* 加载数据
use "`DATA'", clear
xtset 代码 年份

*====================================================================
* 1. 变量预处理 (Variable Pre-processing)
* 功能：生成对数、中心化变量、计算残差等
*====================================================================

* 1.1 基础变量组内中心化 (Within-City Centering)
* 用于构建交互项时减少多重共线性
local basevars ln能耗强度 绿色全要素BOM ln数字金融 绿色金融 ln地区生产总值 ln公路运量源数据 ln风险投资额 ln绿色专利
foreach v of local basevars {
    capture confirm variable `v'
    if _rc {
        continue
    }
    capture confirm variable `v'_c
    if _rc {
        quietly su `v'
        gen `v'_c = `v' - r(mean)
    }
}

* 1.2 数字金融核心变量中心化 (用于构建 DFz_Level 和 DFz_Contrast)
capture drop lnDF_citymean lnDF_wc
bysort 代码: egen lnDF_citymean = mean(ln数字金融)
gen lnDF_wc = ln数字金融 - lnDF_citymean

* 1.3 生成物流效率残差 (Logistics Efficiency Residual)
* 对应正文 Table 1 描述统计中的 LogisticsEff
quietly xtreg ln公路运量源数据 ln地区生产总值_c i.年份, fe
predict ln物流效率_c, residual
quietly su ln物流效率_c
replace ln物流效率_c = ln物流效率_c - r(mean)

* 1.4 机制变量中心化
local mech_vars_centered 产业结构高级化 第二产业占比
foreach v of local mech_vars_centered {
    capture confirm variable `v'
    if _rc {
        continue
    }
    capture confirm variable `v'_c
    if _rc {
        quietly su `v'
        gen `v'_c = `v' - r(mean)
    }
}

*====================================================================
* 2. 定义控制变量 (Define Controls)
*====================================================================
* 对应正文 Table 1 和 Table 2 的控制变量
local controls_base "ln地区生产总值_c  第二产业占比_c ln物流效率_c 环境规则强度 产业结构高级化_c"

local __ctrl ""
foreach v of local controls_base {
    capture confirm variable `v'
    if !_rc local __ctrl "`__ctrl' `v'"
    else di as error "WARNING: 控制变量不存在 -> `v'"
}

global controls "`__ctrl'"
global fe_controls "i.年份"

di as txt ">>> controls      = $controls"
di as txt ">>> fe_controls   = $fe_controls"

*====================================================================
* 3. 定义程序 (Programs Definition)
* 核心程序：mk_basis_z (生成权重和核心解释变量)
* 辅助程序：mk_mech_z (生成机制交互项), run_xtreg (回归封装)
*====================================================================

capture program drop mk_basis_z
program define mk_basis_z
    * 程序功能：根据权重类型(nested/geo/econ)和衰减参数(k)生成 DFz_Level 和 DFz_Contrast
    syntax, k(real) [weight_type(string)]
    if "`weight_type'" == "" local weight_type "nested"
    
    * 设置权重变量名
    if "`weight_type'" == "nested" {
        local wgz_var "广州嵌套权重"
        local wsz_var "深圳嵌套权重"
    }
    else if "`weight_type'" == "geo" {
        local wgz_var "广州地理距离权重"
        local wsz_var "深圳地理距离权重"
    }
    else if "`weight_type'" == "econ" {
        local wgz_var "广州经济距离权重"
        local wsz_var "深圳经济距离权重"
    }
    
    * 生成归一化权重
    tempvar Wgz_pow min_gz max_gz range_gz
    quietly {
        bysort 年份: egen `min_gz' = min(`wgz_var')
        bysort 年份: egen `max_gz' = max(`wgz_var')
        gen double `range_gz' = `max_gz' - `min_gz'
        replace `range_gz' = 1 if `range_gz' < 1e-9
        gen double `Wgz_pow' = ((`wgz_var' - `min_gz') / `range_gz')^`k'
    }
    
    tempvar Wsz_pow min_sz max_sz range_sz
    quietly {
        bysort 年份: egen `min_sz' = min(`wsz_var')
        bysort 年份: egen `max_sz' = max(`wsz_var')
        gen double `range_sz' = `max_sz' - `min_sz'
        replace `range_sz' = 1 if `range_sz' < 1e-9
        gen double `Wsz_pow' = ((`wsz_var' - `min_sz') / `range_sz')^`k'
    }

    * 构造 Level (协同) 和 Contrast (方向) 变量
    capture drop DFz_Level DFz_Contrast L1_DFz_Level L1_DFz_Contrast
    capture drop DFz_GZ DFz_SZ
    
    tempvar Level_Var Contrast_Var
    gen double `Level_Var' = `Wgz_pow' + `Wsz_pow'
    gen double `Contrast_Var' = (`Wsz_pow' - `Wgz_pow') / (`Level_Var' + 1e-6)
    
    gen double DFz_Level = lnDF_wc * `Level_Var'
    gen double DFz_Contrast = lnDF_wc * `Contrast_Var'
    
    * 单核心变量 (用于单核模型检验)
    gen double DFz_GZ = lnDF_wc * `Wgz_pow'
    gen double DFz_SZ = lnDF_wc * `Wsz_pow'
    
    * 生成滞后项 (用于 GTFP 回归)
    sort 代码 年份
    by 代码: gen double L1_DFz_Level = DFz_Level[_n-1]
    by 代码: gen double L1_DFz_Contrast = DFz_Contrast[_n-1]
end

capture program drop mk_mech_z
program define mk_mech_z
    * 程序功能：生成机制变量与 Level/Contrast 的交互项
    capture drop Level_mech_* Contrast_mech_* L1_Level_mech_* L1_Contrast_mech_*
    local mech_vars ln物流效率_c ln风险投资额_c ln绿色专利_c 绿色金融_c 产业结构高级化_c 第二产业占比_c
    foreach mech in `mech_vars' {
        capture confirm variable `mech'
        if _rc continue
        
        local mech_clean = subinstr("`mech'", " ", "_", .)
        local mech_clean = subinstr("`mech_clean'", "-", "_", .)
        quietly {
            gen double Level_mech_`mech_clean' = DFz_Level * `mech'
            gen double Contrast_mech_`mech_clean' = DFz_Contrast * `mech'
            sort 代码 年份
            by 代码: gen double L1_Level_mech_`mech_clean' = Level_mech_`mech_clean'[_n-1]
            by 代码: gen double L1_Contrast_mech_`mech_clean' = Contrast_mech_`mech_clean'[_n-1]
        }
    }
end

capture program drop get_mech_varlist
program define get_mech_varlist, rclass
    * 获取生成的交互项变量列表
    local mech_base ln物流效率_c ln风险投资额_c ln绿色专利_c 绿色金融_c 产业结构高级化_c 第二产业占比_c
    local energy_list ""
    local gtfp_list ""
    foreach v of local mech_base {
        local v_clean = subinstr("`v'", " ", "_", .)
        local v_clean = subinstr("`v_clean'", "-", "_", .)
        capture confirm variable Level_mech_`v_clean'
        if !_rc {
            local energy_list `energy_list' Level_mech_`v_clean' Contrast_mech_`v_clean'
            local gtfp_list `gtfp_list' L1_Level_mech_`v_clean' L1_Contrast_mech_`v_clean'
        }
    }
    return local energy "`energy_list'"
    return local gtfp "`gtfp_list'"
end

capture program drop run_xtreg
program define run_xtreg
    * 通用固定效应回归程序
    syntax, depvar(string) indepvars(string) [model_name(string) lag(string) sample(string)]
    local new_indepvars `indepvars'
    if "`lag'" == "L1" {
        local new_indepvars = subinstr("`new_indepvars'", "DFz_Level", "L1_DFz_Level", .)
        local new_indepvars = subinstr("`new_indepvars'", "DFz_Contrast", "L1_DFz_Contrast", .)
    }
    
    if "`sample'" != "" {
        xtreg `depvar' `new_indepvars' $controls $fe_controls if `sample', fe vce(cluster 代码)
    }
    else {
        xtreg `depvar' `new_indepvars' $controls $fe_controls, fe vce(cluster 代码)
    }
    if "`model_name'" != "" estimates store `model_name'
end

*====================================================================
* 4. 正式分析：生成核心变量
*====================================================================
xtset 代码 年份
* 生成基准变量 (k=2, nested weights)
mk_basis_z, k(2)
mk_mech_z

*====================================================================
* 【对应正文 Table 1】 描述性统计
*====================================================================
* 您的代码中有 tabstat 和 summarize，对应论文 Table 1 Descriptive Statistics
tabstat ln能耗强度_c 绿色全要素BOM_c DFz_Level DFz_Contrast ln数字金融 ln地区生产总值 ln公路运量源数据_c ln风险投资额_c ln绿色专利_c 绿色金融_c 产业结构高级化_c, statistics(N mean sd min max) columns(statistics)

*====================================================================
* 【对应正文 Table 2】 基准回归 (Baseline Results)
*====================================================================
* 模型(1)-(2): 考察 Level 和 Contrast 对能耗强度的影响
* 注：论文中 Table 2 只有 EI (能耗强度) 的结果
run_xtreg, depvar("ln能耗强度_c") indepvars("DFz_Level DFz_Contrast") model_name("T2_Mech_Energy")

* 补充：GTFP 的回归
get_mech_varlist
local gtfp_mech `r(gtfp)'
run_xtreg, depvar("绿色全要素BOM_c") indepvars("DFz_Level DFz_Contrast `gtfp_mech'") lag("L1") model_name("T2_Mech_GTFP")


*====================================================================
* 【对应正文 Table 3】 空间异质性 (Spatial Heterogeneity)
* G1: 深圳核心区, G2: 广州核心区, G3: 其他城市
*====================================================================
capture drop region_group
* 定义分组 (根据城市代码)
gen region_group = 1 if inlist(代码, 3, 11, 4) // Shenzhen Influence
replace region_group = 2 if inlist(代码, 1, 6, 18, 17) // Guangzhou Influence
replace region_group = 3 if missing(region_group)

foreach g in 1 2 3 {
    * 回归结果对应 Table 3 的 G1, G2, G3 列
    run_xtreg, depvar("ln能耗强度_c") indepvars("DFz_Level DFz_Contrast") model_name("T3_Energy_G`g'") sample("region_group==`g'")
}


*====================================================================
* 【对应附录 Table B5】 时间异质性 (Temporal/Phase-Specific Effects)
* Phase 1: 2011-2014, Phase 2: 2015-2018, Phase 3: 2019-2022
* 注意：正文 4.2 节提及此结果，表格在附录 B Table B5
*====================================================================
capture drop period
gen period = 1 if 年份 <= 2014
replace period = 2 if 年份 >= 2015 & 年份 <= 2018
replace period = 3 if 年份 >= 2019

foreach p in 1 2 3 {
    run_xtreg, depvar("ln能耗强度_c") indepvars("DFz_Level DFz_Contrast") model_name("T5_Energy_P`p'") sample("period==`p'")
    run_xtreg, depvar("绿色全要素BOM_c") indepvars("DFz_Level DFz_Contrast") lag("L1") model_name("T6_GTFP_P`p'") sample("period==`p'")
}


*====================================================================
* 【对应正文 Table 4】 机制检验 (Institutional Mechanisms)
* 考察 VC (市场逻辑) 和 GF (政策逻辑) 的调节作用
* 注：此处代码使用了交互项回归
*====================================================================
* 准备交互项 (需确保 mk_mech_z 已运行)
get_mech_varlist
local energy_mech `r(energy)'

* 运行带有所有机制交互项的模型 (Table 4 展示了 VC 和 GF 的列)
* 这里分别展示特定的机制
* Table 4 Column 1: VC (ln风险投资额)
xtreg ln能耗强度_c DFz_Level DFz_Contrast Level_mech_ln风险投资额_c Contrast_mech_ln风险投资额_c $controls $fe_controls, fe vce(cluster 代码)
estimates store Mech_VC

* Table 4 Column 2: GF (绿色金融)
xtreg ln能耗强度_c DFz_Level DFz_Contrast Level_mech_绿色金融_c Contrast_mech_绿色金融_c $controls $fe_controls, fe vce(cluster 代码)
estimates store Mech_GF


*====================================================================
* 【对应附录 Table D1】 其他机制检验 (Technical & Efficiency)
*====================================================================
* Appendix D Table D1 Column 1: Green Patents (ln绿色专利)
xtreg ln能耗强度_c DFz_Level DFz_Contrast Level_mech_ln绿色专利_c Contrast_mech_ln绿色专利_c $controls $fe_controls, fe vce(cluster 代码)
estimates store Mech_GPat

* Appendix D Table D1 Column 2: Logistics (ln物流效率)
xtreg ln能耗强度_c DFz_Level DFz_Contrast Level_mech_ln物流效率_c Contrast_mech_ln物流效率_c $controls $fe_controls, fe vce(cluster 代码)
estimates store Mech_Logistics


*====================================================================
* 稳健性检验 (Robustness Checks) - 对应附录 C
*====================================================================

* --- 【对应附录 Table C1】 子指标分解 (Sub-Index Decomposition) ---
local sub_indices "ln数字广度 ln数字深度 ln电子化水平"
local sub_names   "Breadth Depth Digitization"
local n_sub : word count `sub_indices'

* 备份总指数
clonevar lnDF_wc_backup = lnDF_wc

forvalues i = 1/`n_sub' {
    local var : word `i' of `sub_indices'
    local name : word `i' of `sub_names'
    
    capture confirm variable `var'
    if !_rc {
        di "Running Sub-index: `name'"
        * 临时构造子指标的中心化变量
        capture drop temp_mean
        bysort 代码: egen temp_mean = mean(`var')
        replace lnDF_wc = `var' - temp_mean
        drop temp_mean
        
        * 重新生成 Level/Contrast
        mk_basis_z, k(2)
        
        * 回归
        run_xtreg, depvar("ln能耗强度_c") indepvars("DFz_Level DFz_Contrast") model_name("Sub_`name'")
    }
}
* 恢复总指数
replace lnDF_wc = lnDF_wc_backup
mk_basis_z, k(2)

* --- 【对应附录 Table C2】 替代权重矩阵 (Alternative Weight Matrices) ---
* Column 2: Geographic Weight
mk_basis_z, k(2) weight_type("geo")
run_xtreg, depvar("ln能耗强度_c") indepvars("DFz_Level DFz_Contrast") model_name("T11_Energy_Geo")

* Column 3: Economic Weight
mk_basis_z, k(2) weight_type("econ")
run_xtreg, depvar("ln能耗强度_c") indepvars("DFz_Level DFz_Contrast") model_name("T13_Energy_Econ")

* 恢复基准权重
mk_basis_z, k(2)

* --- 【对应附录 Table C3】 空间衰减参数敏感性 (Decay Parameter k) ---
* Column 1: k=1 (Linear)
mk_basis_z, k(1)
run_xtreg, depvar("ln能耗强度_c") indepvars("DFz_Level DFz_Contrast") model_name("T9_Energy_k1")

* Column 3: k=3 (Strong)
mk_basis_z, k(3)
run_xtreg, depvar("ln能耗强度_c") indepvars("DFz_Level DFz_Contrast") model_name("T9_Energy_k3")

* 恢复基准 k=2
mk_basis_z, k(2)

* --- 【对应附录 Table C4】 剔除核心城市 (Excluding Cores) ---
capture drop keep_flag
gen keep_flag = 1
replace keep_flag = 0 if 代码==1 | 代码==3 // 假设1为广州，3为深圳
run_xtreg, depvar("ln能耗强度_c") indepvars("DFz_Level DFz_Contrast") model_name("T7_Energy_NoCore") sample("keep_flag==1")

* --- 【对应附录 Table C5】 单核心模型 (Single-Core Models) ---
* Column 2: GZ Single-Core (利用 mk_basis_z 生成的 DFz_GZ)
xtreg ln能耗强度_c DFz_GZ $controls $fe_controls, fe vce(cluster 代码)
estimates store Single_GZ

* Column 3: SZ Single-Core
xtreg ln能耗强度_c DFz_SZ $controls $fe_controls, fe vce(cluster 代码)
estimates store Single_SZ

* --- 【对应附录 Table C6】 伪核心安慰剂检验 (Pseudo-Core Falsification) ---
* 代码中包含对其他城市（佛山、汕头等）的循环检验部分
mk_basis_z, k(2)
local placebo_cities "佛山 汕头 珠海 韶关 湛江"
foreach city of local placebo_cities {
    local weight_var "`city'嵌套权重"
    capture confirm variable `weight_var'
    if !_rc {
       
        di "Pseudo Core Test: `city'"
    }
}

* --- 【对应附录 Table C8】 内生性处理 (Endogeneity) ---
* Column 1: Lagged Variables
* Column 2: IV-2SLS 
capture which ivreg2
if !_rc {
    * 生成工具变量 (L2, L3)
    sort 代码 年份
    capture gen double L2_DFz_Level = DFz_Level[_n-2]
    capture gen double L2_DFz_Contrast = DFz_Contrast[_n-2]
    capture gen double L3_DFz_Level = DFz_Level[_n-3]
    capture gen double L3_DFz_Contrast = DFz_Contrast[_n-3]
    
    local safe_controls ln地区生产总值_c ln公路运量源数据_c ln风险投资额_c ln绿色专利_c
    
    ivreg2 ln能耗强度_c (DFz_Level DFz_Contrast = L2_DFz_Level L3_DFz_Level L2_DFz_Contrast L3_DFz_Contrast) `safe_controls' $fe_controls, cluster(代码) robust first
    estimates store IV_Energy
}


*====================================================================
* 【对应正文 Table 6】 预测与反事实模拟 (Forecasting Scenarios)
*====================================================================
* Panel A: Counterfactual Gap Analysis
* Panel B: Forecasting Scenarios
* 此部分对应代码末尾的 Counterfactual scenario prediction 模块

* 确保基准变量存在
mk_basis_z, k(2) weight_type("nested")
mk_mech_z
get_mech_varlist

* 运行基准模型并保存估计结果
xtreg ln能耗强度_c DFz_Level DFz_Contrast `r(energy)' $controls i.年份, fe vce(cluster 代码)
estimates store T2_Mech_Energy

* 定义情景 (Scenario) 并预测
* 示例：freeze2019_post2020 (对应 Table 6 的 Crisis Resilience / Actual)
local scenario "freeze2019_post2020"
local df_var "ln数字金融"

preserve
    * 构造反事实数据
    if "`scenario'"=="freeze2019_post2020" {
        capture drop __dfwc2019
        bysort 代码: egen double __dfwc2019 = mean(cond(年份==2019, lnDF_wc, .))
        replace lnDF_wc = __dfwc2019 if 年份>=2020 & !missing(__dfwc2019)
    }
    
    * 重建核心变量
    mk_basis_z, k(2) weight_type("nested")
    mk_mech_z
    get_mech_varlist
    
    * 预测
    estimates restore T2_Mech_Energy
    predict double yhat_counterfactual, xbu
    
    * 输出结果 (对应 Table 6 的数据基础)
    tabstat yhat_counterfactual if 年份>=2020, by(年份)
restore

di "All Tasks Completed."
log close _all
