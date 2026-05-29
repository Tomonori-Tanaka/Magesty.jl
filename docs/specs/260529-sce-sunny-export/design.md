# Design: SCE → Sunny.jl LSWT エクスポータ

Status: draft (2026-05-29) — M2 スパイクで設計確定

## Summary

フィット済み `SCEModel` を走査して、各 SALC を従来型スピンモデルの相互作用に
分解し、Sunny.jl のフルパイプライン・スクリプト（文字列）を生成する。テキスト
生成のみなので Magesty コアに Sunny 依存は入れない。Sunny が要るのは往復検証
テストのみ（テスト環境）。

**セル戦略（M2 で確定）**: **Magesty 自身の spglib 結果（`Symmetry`）から化学
プリミティブセルを構築**して Sunny `Crystal` に渡す。Sunny は大きなスーパーセル
では対称性解析を無効化する（"Cell is N times too large"）ため、`reshape_supercell`
や Sunny の `primitive_cell` をスーパーセルに対して使う当初案は不可。小さな
プリミティブセルを直接渡せば Sunny は完全な対称性解析を行う（M2 で検証:
febcc→Im-3m, fept→P4/mmm, fege→P2₁3）。Magesty の軌道列挙と同じ spglib 並進を
使うので規約ズレも起きない。

相互作用 `J_ij`（プリミティブセル＋ボンド）は磁気秩序と独立なので、生成スクリプト
は「プリミティブセル＋相互作用」を出力し、**磁気構造（伝播ベクトル `k`／磁気
スーパーセル）はユーザ編集セクション**に分離する。これにより化学 BZ に沿った
折り畳まれない分散が得られ、かつ AFM/スパイラル等の磁気胞にも対応できる（化学
プリミティブ ≠ 磁気単位胞の問題を設計で回避）。

**保証付きフォールバック**: `:explicit`（スーパーセルを P1=全原子別タイプで
`Crystal` 化、`to_inhomogeneous`＋`set_exchange_at!` に最小像オフセット）。M2 で
全 fixture を `energy(sys)` と `predict_energy-j0` で機械精度一致を確認済み。
ただし磁気胞=スーパーセルで分散は折り畳まれる。`placement = :primitive`（既定）/
`:explicit` で切替。

### エネルギー規約（M2 実測で確定）

- スピン `s=1`、`Moment(s=1, g=2)`、`:dipole` モード。
- 交換: `set_exchange*`/`Bond` に 3×3 行列 `M` を**そのまま**渡す
  （`energy(sys)` の交換項 = `Σ ê_a' M ê_b`）。向きは `(a,b)`、`a<b` 正準。
- 単イオン: 演算子 `op = c · Σ_{αβ} A[α,β] Sα Sβ`、係数
  `c = 2/(s(2s-1))`（s=1 → 2）。spin 系の二次演算子の古典期待値が
  `s(2s-1)/2 · ê'Aê` になるための補正（M2 で s=1,1.5,2 で確認）。
- `j0` と定数 SALC は Sunny に定数項がないため落とす（分散には無影響）。

### プリミティブ展開の正確な成立条件（M2 で実証）

Magesty は **cutoff < L/2**（それより遠い相互作用は無視）を前提に fit する。この
前提下ではプリミティブ展開は厳密に定義でき、reshape で学習スーパーセルへ戻すと
`predict_energy - j0` を機械精度で再現する（M2: dimer と **chain**（ntran=2,
cutoff=L/2）で 1e-16）。

**厳密成立の判定基準**: 全 pair で `multiplicity × #clusters == ntran`。
- cutoff = L/2（共線 ±像の折り重なり）→ `mult=2, #clusters=ntran/2`,
  product = ntran → **クリーン**（問題なし）。
- `product > ntran`（cutoff > L/2 で非共線の異方ボンドが同一原子対に重なる）→
  異なる方向＝異なるテンソルの和になり分離不能 → **真に lossy**。この場合のみ
  ガードで警告し、`:explicit`（折り畳み）へフォールバックする。

**実装の要点（M2 で確定した正しい扱い）**:
- プリミティブボンドの結合は `jphi · scaling · M0`（**multiplicity を掛けない**）。
  ±像（multiplicity）の複製は Sunny の bond-reverse と reshape が自動再現する。
  multiplicity を畳み込むと二重計上になる。
- ボンドオフセットは**最小像**（周期方向のみ wrap）で計算。
- 対称性伝播を避けるため `Crystal(latvecs, positions, 1; ...)` で **P1 を強制**
  （スーパーセルでは Sunny が対称性解析を無効化するため、Magesty の spglib から
  得たプリミティブセルを渡す。M2 検証: febcc→Im-3m, fept→P4/mmm, fege→P2₁3 と
  Sunny は正しく同定するが、相互作用配置は P1 強制＋明示配置が確実）。
- DM（奇 Lf）が境界自己ペアに乗る成分は Magesty が `filter_basisdict` で既に
  除去済み（±像で打ち消すため）。偽の DM は出ない。

検証は `reshape_supercell` で学習スーパーセルへ戻し、ランダム配置で
`energy(sys)` と `predict_energy - j0` を突合（往復テストの門番）。

## Module layout

| Target | Change |
|---|---|
| `src/SunnyExport.jl` (new) | `VaspConvert.jl` と同様にモジュールにせず直接 include。`sce_to_sunny` と内部分解 `_sunny_build_terms` ほか `_sunny_*` ヘルパ。Sunny 非依存。 |
| `src/Magesty.jl` | `include("SunnyExport.jl")`（`VaspConvert.jl` の直後）、`export sce_to_sunny`。 |
| `cli/src/MagestyCLI.jl` | `@cast module Sunny` 追加（`script` サブコマンド）。 |
| `test/component/test_sunny_export.jl` (new) | Sunny 非依存の単体テスト。 |
| `test/integration/sunny_export/test_roundtrip.jl` (new) | Sunny ゲート付き往復テスト。 |
| `docs/src/sunny_export.md` (new) | 使い方＋変換表。 |

## API

```julia
function sce_to_sunny(
    model::SCEModel;
    output::Union{AbstractString,Nothing} = nothing,
    symprec::Real = 1e-5,
    placement::Symbol = :homogeneous,        # :homogeneous | :explicit
    qpath_labels::Union{Nothing,AbstractVector} = nothing,  # future hook
)::String
```

スクリプト全文の文字列を返し、`output` 指定時は `.jl` を付与して書き出す
（`vasp_to_extxyz` と同様）。標準 Julia docstring（# Arguments / # Returns /
# Examples）を付ける。

内部中間表現（Sunny 型を使わない素の struct）:

- `PairTerm`: `i, j, offset::NTuple{3,Int}, J::SMatrix{3,3,Float64}`
- `OnsiteTerm`: `i, A::SMatrix{3,3,Float64}`（対称・トレースレス、`E = êᵀAê`）
- `SkippedTerm`: 警告用の説明文字列

## Types and conventions

分解は `predict_energy` を厳密再現する。オラクルは
`Fitting._predict_energy`（`src/Fitting.jl:1059-1085`）と
`design_matrix_energy_element`（`src/Fitting.jl:762-817`）。

`model.basis.salcbasis.salc_list` を index `ν`（`model.jphi[ν]` と一致）で走査:

1. `n_C = length(group[1].atoms)`、`scaling = Fitting._cluster_scaling(n_C)`
   を**再利用**（`(4π)^(n_C/2)` を再導出しない）。
2. `group[1].ls` / `Lf` で分類:
   - `ls == [1,1]`（2体）→ bilinear pair。採用。
   - `ls == [2]`（1体・`Lf==2`）→ single-ion。採用（`Lf=0` は定数で `j0` に
     吸収、奇数 `l` はパリティ禁制）。
   - それ以外 → `SkippedTerm`。
3. **3×3 を `folded_tensor` から直接構築**（式非依存で厳密）。`l1=l2=1` では
   `M[μ,ν] = (3/4π)·permute(folded_tensor)`、テッセラル添字写像 `m_idx 1→y,
   2→z, 3→x`（`Z_{1,m}=√(3/4π)·e_μ` 由来）。これはカーネルが行う縮約から `Z` の
   前因子を括り出したものに一致し、項ごとに厳密再現する。`technical_notes.md` の
   `J=√3·J_0` / `D=(3/√2)·J_1` / `Γ(...)` 式は単体テストで突き合わせ
   （式は**文書上の解釈**、`folded_tensor` 経路が**実装の門番**）。
4. ペア集約は**2段ループ**で行う（`design_matrix_energy_element`
   `src/Fitting.jl:786-814` に対応）。1 つの key_group には軌道内のプリミティブ
   原子ごとに**複数の `cbc`** が入り、各 `cbc` は自身の `multiplicity`（スカラの
   重み）と `clusters`（並進像のリスト）を持つ。SALC ごとに 1 ボンドだけ取る実装
   にすると残りのプリミティブ原子を取りこぼすため、必ず両ループを回す:

   ```
   for cbc in key_group              # 軌道内プリミティブ原子ごと（複数）
       for cluster in cbc.clusters   # 並進像ごと（複数のボンド）
           a, b = cluster[1], cluster[2]            # site 順（sorted ではない）
           sa, sb = subl(a), subl(b)                # map_s2p で副格子 index
           n = offset(a, b)                         # 下記「ボンドオフセット」参照
           # Sunny 正準ボンドは sublattice 昇順。向きが一致すれば M、
           # 反転していれば Mᵀ（Lf=1 の DM は反転で符号反転するため必須）。
           if (sa, n) <= (sb, ...)                  # 正準向き判定
               accumulate!(dict, Bond(sa, sb, n),  jphi[ν]*scaling*cbc.multiplicity*M)
           else
               accumulate!(dict, Bond(sb, sa, -n), jphi[ν]*scaling*cbc.multiplicity*transpose(M))
           end
       end
   end
   ```

   `M` は cluster の `(a,b)` に対する 3×3（手順 3）。`enumerate_orbit_clusters` は
   sorted multiset で重複排除するため無順序ペアは軌道内で一度だけ＝内部二重計上は
   ないが、**site 順は sorted ではない**ので、並進後に `a` の原子 index が `b` を
   上回ることがある。正準向き判定と転置はこのために必要。異なる key_group が同一
   正準ボンドに落ちる場合は `Dict` 上で合算される。
5. **単イオン**: rank-1 `folded_tensor`（`m=-2..2`）を `Z_{2,m}` の Cartesian
   多項式で縮約 → 対称トレースレス `A`。`jphi[ν]·scaling·multiplicity` 倍して
   副格子ごとに加算。`set_onsite_coupling!` に `spin_matrices` 由来のスピン演算
   多項式として出力（`S` 正規化の前因子は往復テストで固定）。

### ボンドオフセット / 副格子写像

- `:explicit`: サイト＝スーパーセル原子の絶対 index、offset 常に `(0,0,0)`、
  `set_exchange_at!` 使用。オフセット計算不要。
- `:homogeneous`: スーパーセル上に `Crystal` を作る。副格子 index ＝原子 index
  （Sunny が並べ替える場合に備え位置→サイトの逆引きを 1 度だけ構築）。物理的な
  ボンド変位（整数 `offset`）は位置から最小像
  `supercell.lattice_vectors * (x_frac[:,b]-x_frac[:,a]+n)` で復元、`offset=n`。
  その後 `reshape_supercell(sys, primitive_cell(cryst))` でプリミティブへ縮約。
  `Symmetry.map_s2p` / `map_p2s` / `atoms_in_prim` は照合用に留める。

### 生成スクリプトの構成

1. ヘッダコメント: 出自、`s=1, g=2`、単位 eV、**スキップ項一覧**。
2. `using Sunny, GLMakie`。
3. `Crystal(latvecs, positions; types, symprec)`。
4. `System(cryst, [i => Moment(s=1, g=2) ...], :dipole; dims=(1,1,1))`。
5. 相互作用: `set_exchange!`/`Bond`（homogeneous）or `set_exchange_at!`
   （explicit）、`set_onsite_coupling!`。行列はリテラル `[a b c; d e f; g h i]`。
6. `sys2 = reshape_supercell(sys, primitive_cell(cryst))`（homogeneous のみ）＋
   `minimize_energy!(sys2)`（`randomize_spins!` リトライ例はコメント）。
7. `SpinWaveTheory(sys2; measure=ssf_perp(sys2))`、**`# EDIT ME`** の q-path
   ブロック（`qs`/`labels` 例）、`q_space_path`、`dispersion`。
8. プロット（`plot_intensities` / `lines`）。

## Impact on linked sites

- [ ] Spherical-harmonics convention (`TesseralHarmonics`): 読み取りのみ（`Z_{1,m}`
      / `Z_{2,m}` の Cartesian 形を変換に使用）。規約変更なし。
- ~~SCE coefficient XML (`save` / `load`)~~: 影響なし（読み書きとも変更しない）。
- [ ] `Fitting` ↔ `SALCBasis`: 読み取りのみ（`salc_list` の順序・`folded_tensor` /
      `multiplicity` / `clusters` を使用）。
- [ ] `.claude/agents/` references: モジュール追加に伴い必要なら test 対象を追記。
- [ ] `SPEC.md` / `docs/src/api.md` updates: 新公開関数を反映。`docs/make.jl` の
      `pages` には Tools セクション（`tools.md` と同列）に `sunny_export.md` を追加。

## Test strategy

- **単体（Sunny 非依存）** `test/component/test_sunny_export.jl`:
  Heisenberg-only → `M==J·I`、DM-only → 反対称、Γ-only → 対称トレースレスを
  `technical_notes.md` の値と照合。スキップ判定（3体 / `l=3`）。ボンド正準化・
  マージ（転置）。単イオン `A` の対称性・`tr(A)=0`。文字列の節存在と
  `Meta.parseall` での構文妥当性。
- **往復（Sunny ゲート付き）**: 既存の integration レイアウト（fixture ごとの
  サブディレクトリ）に合わせて配置する（例 `test/integration/fept_tetragonal_2x2x2/`、
  `test/integration/dimer/`）。`Sunny` をテスト環境のみに追加（コア `Project.toml`
  に入れない）、`MAGESTY_TEST_SUNNY` などで CI ゲート。`_build_terms` の IR から
  System を直接再構築し、≥20 個のランダム単位スピンで
  `energy(sys) ≈ predict_energy(model) - j0` を確認。fixture は dimer /
  fept_tetragonal_2x2x2 / 異方 fixture（fege 系があれば）。ここで `placement`
  既定を確定し、`reshape_supercell` の homogeneous 動作も確認。

## Risks and open items

1. **Sunny の空間群伝播 vs Magesty の並進軌道テンソル**（異方 J/DM）: `:homogeneous`
   の最重要リスク。往復テストで判定、`:explicit` を保証付きフォールバックに。
2. **単イオンの `S` 正規化**（`spin_matrices(1)` の古典還元）: 往復テストで数値固定。
3. **inhomogeneous 系の `reshape_supercell`**: 非対応と想定。スパイクで確認。
4. **Sunny の原子並べ替え**: 位置→サイト逆引きでガード。
5. **`j0` のドロップ**: Sunny に定数項なし。絶対エネルギーは `j0` だけずれる
   （分散には無影響）。ドキュメント明記。
6. **DM の符号規約**: 生の 3×3（`+` 規約、Magesty と一致）を出力。`-J` 文献規約は
   ユーザ側で反転する旨を明記。

### 実装順（tasklist と対応）

spec 合意 → サブモジュール骨組み → 早期 Sunny スパイク → 分解 → オフセット＋
文字列生成 → wire-in/export → 単体テスト → Sunny テスト環境＋往復テスト →
CLI → docs。
