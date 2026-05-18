# Regression test for the energy-centered-design-matrix refactor.
#
# The reference `(j0, jphi)` literals below were captured from the
# pre-refactor bias-in-`X` code (see
# `tools/personal/capture_energy_centered_refs.jl`,
# `tools/personal/capture_fege_ridge_w10.jl`). They pin the
# pre-refactor floating-point output for representative
# `(estimator, weight)` pairs on two integration datasets.
#
# Case set:
# - `febcc_2x2x2_pm`: well-conditioned (`num_spinconfigs > num_salcs`),
#   so all four `(estimator, weight)` pairs
#   `{OLS, w in {0.0, 0.5, 1.0}}` and `Ridge(lambda = 0.1) at w = 0.5`
#   give numerically unique coefficient vectors.
# - `fege_2x2x2`: `num_salcs >= num_spinconfigs`, so the energy block
#   alone is rank-deficient. Energy-weighted fits (`weight < 1`) admit
#   a non-unique minimum-norm solution whose specific point in the
#   null space depends on how the QR factorization breaks ties; the
#   pre-refactor (bias-in-`X`) and post-refactor (centered) systems
#   pick different but equally valid points and the corresponding
#   `(j0, jphi)` numbers legitimately disagree at O(1). Only
#   torque-weight-one cases are retained here, where the torque block
#   resolves the rank deficiency and the solution is unique up to
#   floating-point rounding noise.
#
# The refactor (energy-only centering inside
# `assemble_weighted_problem` + `extract_j0_jphi` with no bias-column
# discard) is mathematically equivalent to the bias-in-`X`
# formulation, but the floating-point operation sequence differs.
# Agreement is therefore expected within rounding noise rather than
# bit-for-bit. `REGRESSION_RTOL` is measured below.

using Test
using TOML
using Magesty

# Regression tolerance. Set from the maximum observed `rtol` across
# all retained cases, padded by 10x and rounded up to the next power
# of 10. The pre-refactor (bias-in-`X`) and post-refactor (centered)
# formulations are mathematically equivalent for these cases, so
# observed differences are pure floating-point rounding noise from
# the different operation sequences (LAPACK QR on a matrix of
# different shape; `mean(...)` summation order; etc.).
#
# Measured maxima (this commit):
# | dataset           | case        | rel(j0)  | rel(jphi)max  |
# |-------------------|-------------|----------|---------------|
# | febcc_2x2x2_pm    | OLS_W00     | 4.4e-16  | 3.8e-12       |
# | febcc_2x2x2_pm    | OLS_W05     | 6.6e-16  | 2.3e-11  <-- worst
# | febcc_2x2x2_pm    | OLS_W10     | 0.0      | 1.3e-15       |
# | febcc_2x2x2_pm    | RIDGE01_W05 | 2.2e-16  | 6.9e-12       |
# | fege_2x2x2        | OLS_W10     | 0.0      | 9.5e-13       |
# | fege_2x2x2        | RIDGE01_W10 | 0.0      | 2.4e-13       |
#
# Worst observed = 2.3e-11; padded 10x = 2.3e-10; rounded up to the
# next power of 10 = 1e-9.
const REGRESSION_RTOL = 1e-9

# --- Reference values (pre-refactor, captured M1) -------------------

# === FEBCC2X2X2PM_OLS_W00 ===
# j0:
const J0_FEBCC2X2X2PM_OLS_W00 = -128.66348641968855
# jphi (length = 7):
const JPHI_FEBCC2X2X2PM_OLS_W00 = Float64[
    0.003214331407123817,
    0.004874446323422749,
    -0.0004626704452941751,
    0.0006149062347903999,
    -0.0010469467470516137,
    -0.06710689988233057,
    -0.002510081212154956,
]

# === FEBCC2X2X2PM_OLS_W05 ===
# j0:
const J0_FEBCC2X2X2PM_OLS_W05 = -128.66298661864556
# jphi (length = 7):
const JPHI_FEBCC2X2X2PM_OLS_W05 = Float64[
    0.003136200699199452,
    0.004332277787658181,
    -0.00044534497920633876,
    0.0007582874246008081,
    -0.0009535333919443681,
    -0.06753550360583406,
    -0.0025342969805366962,
]

# === FEBCC2X2X2PM_OLS_W10 ===
# j0:
const J0_FEBCC2X2X2PM_OLS_W10 = -128.65031500146242
# jphi (length = 7):
const JPHI_FEBCC2X2X2PM_OLS_W10 = Float64[
    0.0028103502442177065,
    0.0011994027578327202,
    0.00035232011958135923,
    0.0010055600751340812,
    -0.00028500299286734717,
    -0.06717280486466949,
    0.0006456721353788952,
]

# === FEBCC2X2X2PM_RIDGE01_W05 ===
# j0:
const J0_FEBCC2X2X2PM_RIDGE01_W05 = -128.66163732492723
# jphi (length = 7):
const JPHI_FEBCC2X2X2PM_RIDGE01_W05 = Float64[
    0.003289415852314114,
    0.004849840227709082,
    -0.00041778094268365203,
    0.0005454935019589928,
    -0.0008894612202714528,
    -0.0656724778491685,
    -0.0019111038493075754,
]

# === FEGE2X2X2_OLS_W10 ===
# j0:
const J0_FEGE2X2X2_OLS_W10 = -412.42766841427454
# jphi (length = 146):
const JPHI_FEGE2X2X2_OLS_W10 = Float64[
    -0.00029196023618591725,
    -1.887321325390358e-5,
    4.694172114086719e-5,
    5.661696525280528e-5,
    -5.941769819151026e-6,
    -3.287640238111194e-5,
    1.914295542511932e-5,
    0.0004329128800810286,
    -0.00016798228365466438,
    0.00042916161110312433,
    9.41332073772409e-5,
    -0.00038594911141989416,
    0.000129417548296919,
    -0.0002555678326711127,
    1.4403502746700679e-5,
    -0.00019741489940243816,
    0.0002117629328676287,
    0.0004662505575145001,
    -0.0001487585486140732,
    -0.00029353897042016207,
    8.54667118629624e-5,
    5.224242041580281e-5,
    0.0001254253617544558,
    3.961491838236033e-5,
    0.0023471213509719774,
    -0.00024417594567059696,
    1.4324872842957706e-5,
    -0.0004511291889100065,
    6.486675226291513e-5,
    0.0002726356543519504,
    -0.00023539668299567788,
    -0.00010005616958868314,
    -0.00016266094939804913,
    0.00012367546244653786,
    -7.38248511811667e-5,
    -9.0679900761417e-5,
    0.0070229687213375054,
    -7.582227232819429e-5,
    0.0011618296338394725,
    -8.978980349974243e-5,
    -9.224133862896864e-5,
    0.00012199167743459809,
    -0.0001034210709416387,
    -1.9536762352313293e-5,
    -0.00023251710167901255,
    -0.00015747064463293554,
    -5.043212247190838e-5,
    -1.4775344926559227e-5,
    -0.00015186103405873024,
    0.00022483963887756091,
    0.00035991997384879034,
    -0.00012416821042800585,
    -0.00022498421654553767,
    6.864799961965702e-5,
    0.005124848821932029,
    -0.00022706207751317686,
    8.515104776327399e-6,
    0.000596485877430042,
    3.701328894144155e-5,
    0.0001971194735458476,
    0.0003869595952054986,
    -9.904624624122809e-5,
    -8.340291687124052e-5,
    0.00040770475615773023,
    0.0001428533774832835,
    -4.0362173146245616e-5,
    0.00016387443714056054,
    0.00012988465641022712,
    0.0003192813820932172,
    0.0001310812419755022,
    -9.32528897390375e-5,
    -0.0001902862833303923,
    0.00388038263054651,
    -0.00024821443394254934,
    -0.00012826084049789906,
    -6.834270798378396e-5,
    6.176541568625006e-5,
    -0.000131341127510928,
    0.00032733213892809527,
    -0.00019369949037989463,
    0.00020117764976412784,
    -8.265504958857486e-5,
    -6.360806635456318e-5,
    -3.074886283639519e-5,
    -6.796672205496087e-5,
    -7.346222804860978e-5,
    -8.784560329500458e-5,
    -0.00026051593113804385,
    -0.0003498109256163143,
    5.3309824477423e-5,
    0.00022615627309974058,
    0.00011363299899240167,
    -0.00013949224264420742,
    -5.443524851288777e-5,
    0.00017346267535562802,
    5.25506058533811e-5,
    0.00020907637085496465,
    -6.975410922909507e-5,
    0.00020486863124333682,
    -0.00035528074947278914,
    1.5722184182814924e-6,
    -0.03828427975139577,
    -0.0007211628556422602,
    0.0007091021610707674,
    -0.0001422676733464037,
    -7.104025197430776e-5,
    -0.00030919868762416106,
    -5.6142225305206667e-5,
    -0.00024446620010041507,
    -0.00019597181636498505,
    0.0020213193551076292,
    -0.00018658130015248596,
    -0.00019905690876497874,
    -1.7250765014544513e-5,
    0.00011291598874350813,
    0.00019222791009942085,
    3.3433307377947365e-5,
    0.00012297623411803,
    0.00020265689870713622,
    -3.9716342320139016e-5,
    -5.980259258355868e-5,
    -4.76242997204596e-6,
    -0.0002424077736808621,
    -1.2366528574100342e-5,
    0.00019093289962188935,
    -0.00021754745748597606,
    -1.0383744057400388e-5,
    6.82673238274905e-5,
    -0.006361243630709884,
    -5.7308216976533434e-5,
    0.00011826774222262535,
    -0.0001953258057194098,
    -0.00012346033377283575,
    -0.00016683444131633776,
    3.8013091989883806e-5,
    0.0002956687987126168,
    2.3573176241910334e-5,
    -0.0002472268220623351,
    3.721064928170928e-5,
    -0.00011926458901779152,
    0.0002399067741609267,
    7.668196972077184e-5,
    7.334734117175541e-7,
    0.0002696941644642969,
    9.067012633954311e-5,
    -6.210279767493177e-5,
]

# === FEGE2X2X2_RIDGE01_W10 ===
# j0:
const J0_FEGE2X2X2_RIDGE01_W10 = -412.4749115051957
# jphi (length = 146):
const JPHI_FEGE2X2X2_RIDGE01_W10 = Float64[
    -0.0018647912902778102,
    0.00013766441703777155,
    -4.5454778207693756e-5,
    5.55295743835704e-5,
    -1.94833952554547e-5,
    -0.00023354433243450814,
    0.0006135768303082948,
    7.864789652434846e-5,
    -0.00011232586591115688,
    9.500917980420792e-5,
    -1.9597111357128033e-5,
    -5.593958266908328e-5,
    6.651750690501165e-5,
    -4.8620354324990096e-5,
    6.651425305246183e-5,
    0.0005813029481625321,
    8.090322760351064e-5,
    5.2305016617514685e-5,
    -3.46162624151687e-5,
    -0.00016244448406408393,
    -1.0636599176450598e-5,
    -7.495503557806108e-5,
    0.0001246098887306896,
    4.303415491234575e-5,
    0.0010558049302765097,
    -6.588520814966952e-5,
    3.6857991121223114e-5,
    0.00047843857328798504,
    6.496260121405798e-5,
    2.0003686922242796e-5,
    3.937419780432861e-6,
    1.1569335051152962e-5,
    -0.00012519924135254477,
    7.898961560922598e-5,
    4.618724235227247e-5,
    0.00018423897611704825,
    -0.0008510081167000985,
    -7.588356713717201e-5,
    0.0001400640643028276,
    -1.038816966665543e-5,
    -0.00019197340136750904,
    3.076304906212954e-5,
    -2.5563683161476057e-6,
    -0.00010011474623102098,
    -0.00010785313599909066,
    0.000696284987309455,
    -8.292182953062725e-5,
    -7.546668144478464e-6,
    -8.871429762608286e-5,
    -2.2792063137410706e-5,
    0.00020679402234318424,
    -3.109555196498709e-5,
    -0.00019286547150046748,
    -3.061523316967502e-5,
    -0.0014309338328533806,
    -0.00013021231987229938,
    -0.00011028772866726113,
    0.00011427639045311349,
    -6.758778226197936e-5,
    -0.00011482984503069263,
    -5.978864455392365e-5,
    -0.00018322052832332852,
    -0.00025471848483216066,
    0.0007536096104153604,
    6.276387543022657e-6,
    -5.8149326081699194e-5,
    0.00017847988629510162,
    -2.5799991737567974e-6,
    2.2326130422254145e-5,
    -5.2335127688588584e-5,
    9.575643719918679e-6,
    7.707461737783797e-5,
    0.0018635058285043826,
    -9.295234296200736e-5,
    -5.891731199553172e-5,
    -5.407993895609365e-6,
    -4.0609859829036e-5,
    -7.135501140125937e-5,
    0.00017626322051906196,
    -0.00019146544654800746,
    5.7694428568721656e-5,
    0.0006195844867214481,
    -3.0273368244063653e-5,
    3.859440955664252e-5,
    -9.527243616298463e-5,
    -9.963562443582575e-5,
    7.019483818683606e-5,
    -4.883357529421767e-5,
    -0.0002017618376501397,
    3.0166750413246698e-5,
    -0.0023925066423314077,
    -0.00010780870620510651,
    -1.5212150948210189e-5,
    -4.6517261383000736e-5,
    0.00010924510503218926,
    1.8890896002949716e-5,
    6.645207046381063e-5,
    -0.00019072864849713312,
    0.0001871682731516128,
    -0.0009564394344558292,
    -0.00010445624824145786,
    -0.01472106780508398,
    -0.00013166832116599708,
    1.2260537723931557e-6,
    -0.00018574951067139645,
    -7.658355638162108e-5,
    0.0002324213608869067,
    -1.8250159433577187e-5,
    -0.00010696160297227875,
    -1.0703529313162996e-5,
    0.001217021808437698,
    -2.1565493293892155e-5,
    2.2578061270856223e-5,
    -0.00014013589123298658,
    -2.7299405972166837e-5,
    2.0753145227842104e-5,
    -4.6320310512555423e-5,
    -8.950156376668275e-6,
    0.00010922232725316209,
    0.0005947372433989102,
    -9.296694867432088e-6,
    2.7165832797433375e-5,
    -1.7744983500419175e-5,
    7.805762164276469e-5,
    8.87008054507222e-5,
    -0.0002506149665963986,
    -4.428307569575231e-5,
    5.457296634256216e-5,
    -0.005186781407848798,
    -6.903159901213483e-5,
    0.00023951968997795844,
    2.348383370330813e-5,
    8.837454843662426e-7,
    -0.0002531508096458943,
    0.0006492450226840426,
    0.0001164575578930001,
    3.9506137502483227e-5,
    -9.325877731242671e-5,
    -6.617233733721958e-5,
    -9.98109239952503e-5,
    8.362670402562625e-5,
    -7.283776272431011e-6,
    -9.72793213420831e-5,
    0.000499677863231159,
    4.5289052516070465e-5,
    -2.7154455581790976e-5,
]

# --- Dataset / fit helpers ------------------------------------------

function _ec_build_dataset(tag::String)
    input_path = joinpath(@__DIR__, "..", "integration", tag, "input.toml")
    embset_path = joinpath(@__DIR__, "..", "integration", tag, "EMBSET")
    input = TOML.parse(open(input_path, "r"))
    basis = SCEBasis(input; verbosity = false)
    spinconfigs_all = Magesty.read_embset(embset_path)
    num_cfg_target = 2 * length(basis.salcbasis.salc_list)
    num_cfg = min(length(spinconfigs_all), num_cfg_target)
    return SCEDataset(basis, spinconfigs_all[1:num_cfg])
end

function _ec_fit(dataset, estimator, w)
    return fit(SCEFit, dataset, estimator; torque_weight = w, verbosity = false)
end

# --- Test cases -----------------------------------------------------

# Each entry: (label, estimator_thunk, weight, j0_ref, jphi_ref).
# The estimator is wrapped in a thunk so the same case tuple format
# can be reused if assertions are ever wrapped in `@test_skip` again.
const _FEBCC_CASES = (
    ("OLS_W00",     () -> OLS(),                  0.0, J0_FEBCC2X2X2PM_OLS_W00,     JPHI_FEBCC2X2X2PM_OLS_W00),
    ("OLS_W05",     () -> OLS(),                  0.5, J0_FEBCC2X2X2PM_OLS_W05,     JPHI_FEBCC2X2X2PM_OLS_W05),
    ("OLS_W10",     () -> OLS(),                  1.0, J0_FEBCC2X2X2PM_OLS_W10,     JPHI_FEBCC2X2X2PM_OLS_W10),
    ("RIDGE01_W05", () -> Ridge(; lambda = 0.1),  0.5, J0_FEBCC2X2X2PM_RIDGE01_W05, JPHI_FEBCC2X2X2PM_RIDGE01_W05),
)

const _FEGE_CASES = (
    ("OLS_W10",     () -> OLS(),                  1.0, J0_FEGE2X2X2_OLS_W10,     JPHI_FEGE2X2X2_OLS_W10),
    ("RIDGE01_W10", () -> Ridge(; lambda = 0.1),  1.0, J0_FEGE2X2X2_RIDGE01_W10, JPHI_FEGE2X2X2_RIDGE01_W10),
)

@testset "Energy-centered design matrix regression" begin
    @testset "febcc_2x2x2_pm" begin
        dataset = _ec_build_dataset("febcc_2x2x2_pm")
        for (label, est_thunk, w, j0_ref, jphi_ref) in _FEBCC_CASES
            f = _ec_fit(dataset, est_thunk(), w)
            @test isapprox(intercept(f), j0_ref; rtol = REGRESSION_RTOL)
            @test isapprox(coef(f), jphi_ref; rtol = REGRESSION_RTOL)
        end
    end

    @testset "fege_2x2x2" begin
        dataset = _ec_build_dataset("fege_2x2x2")
        for (label, est_thunk, w, j0_ref, jphi_ref) in _FEGE_CASES
            f = _ec_fit(dataset, est_thunk(), w)
            @test isapprox(intercept(f), j0_ref; rtol = REGRESSION_RTOL)
            @test isapprox(coef(f), jphi_ref; rtol = REGRESSION_RTOL)
        end
    end
end
