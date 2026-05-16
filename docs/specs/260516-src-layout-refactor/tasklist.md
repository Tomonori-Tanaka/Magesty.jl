# Tasklist: `src/` layout refactor

進行管理は in-session の `TaskCreate` でも併用する。本ファイルは
コミット単位のマイルストーンを記録する。

## Milestone 1: ファイル名統一 + ディレクトリ平坦化

- [ ] `src/common/version.jl` → `src/Version.jl`（`git mv`）
- [ ] `src/common/SortedCounter.jl` → `src/SortedCounters.jl`
- [ ] `src/types/AtomCells.jl` → `src/AtomCells.jl`
- [ ] `src/types/Basis.jl` → `src/Basis.jl`（このコミットでは改称しない）
- [ ] `src/utils/SphericalHarmonicsTransforms.jl` → `src/SphericalHarmonicsTransforms.jl`
- [ ] `src/utils/AngularMomentumCoupling.jl` → `src/AngularMomentumCoupling.jl`
- [ ] `src/utils/ConfigParser.jl` → `src/ConfigParser.jl`
- [ ] `src/utils/atomsbase_adapter.jl` → `src/AtomsBaseAdapter.jl`
- [ ] `src/utils/RotationMatrix.jl` → `src/RotationMatrix.jl`
- [ ] `src/utils/MySphericalHarmonics.jl` → `src/MySphericalHarmonics.jl`（このコミットでは改称しない）
- [ ] `src/utils/xml_io.jl` → `src/XMLIO.jl`
- [ ] 空になった `common/`, `types/`, `utils/` を削除
- [ ] `src/Magesty.jl` の `include(...)` パスを全て top-level に更新
- [ ] `make test-all` / `make test-tools` / `make test-aqua` / `make build` 緑確認
- [ ] commit: `refactor(src): flatten directory structure and unify file naming`

## Milestone 2: モジュール改称

- [ ] `MySphericalHarmonics.jl` → `TesseralHarmonics.jl` + `module` 名変更
- [ ] `Basis.jl` → `CoupledBases.jl` + `module` 名変更
- [ ] `Optimize.jl` → `Fitting.jl` + `module` 名変更
- [ ] src 内の `using ..Basis` → `using ..CoupledBases` 全置換
- [ ] src 内の `Basis.CoupledBasis` 等のドット参照を `CoupledBases.CoupledBasis` に置換
- [ ] src 内の `MySphericalHarmonics` 参照を `TesseralHarmonics` に置換
- [ ] src 内の `Optimize` モジュール参照を `Fitting` に置換
- [ ] `src/Magesty.jl` の `include` / `using .X` 行を新名で更新
- [ ] test の `using` 文 + ファイル名（`test_MySphericalHarmonics.jl` →
      `test_TesseralHarmonics.jl`）更新
- [ ] tools/ の `include` / `using` の追随修正（grep で確認）
- [ ] `docs/src/api.md` の `@docs` ブロック内モジュール名を新名に
- [ ] `make test-all` / `make test-tools` / `make test-aqua` /
      `make test-jet` / `make build` 緑確認
- [ ] commit: `refactor(src): rename internal modules for clarity`

## Milestone 3: ドキュメント同期

- [ ] `SPEC.md` のディレクトリ図と「主要モジュール」「ユーティリティ」
      表を新レイアウト・新名前で書き直し
- [ ] `CLAUDE.md` 物理規約セクション「球面調和関数の規約」内、
      `MySphericalHarmonics` 名指しを `TesseralHarmonics` に
- [ ] `docs/design-notes/` 内の旧名参照を最小限で同期
- [ ] `git grep` で以下がゼロ:
      - `MySphericalHarmonics`
      - `module Basis\b`
      - `module Optimize\b`
      - `common/version`, `common/SortedCounter`
      - `types/AtomCells`, `types/Basis`
      - `utils/atomsbase_adapter`, `utils/xml_io`,
        `utils/MySphericalHarmonics`, `utils/SphericalHarmonicsTransforms`,
        `utils/AngularMomentumCoupling`, `utils/ConfigParser`,
        `utils/RotationMatrix`
- [ ] commit: `docs: sync SPEC.md and CLAUDE.md with src layout refactor`

## Milestone 4: PR

- [ ] PR description に「公開 API 不変」と「改称マッピング表」を明記
- [ ] reviewer に公開準備の文脈を伝達

## 検証コマンド（各 milestone 後に共通実行）

```bash
make test-all
make test-tools
make test-aqua
make test-jet
make build
```
