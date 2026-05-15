# Config4System の役割の見直し

**Status**: 未着手（2026-05-14 提案）
**関連**: [`docs/specs/260514-sce-public-api/`](../specs/260514-sce-public-api/)

## 背景

SCE 公開 API リファクタ（4 型構成）の AtomsBase 統合を実装した際、
adapter (`src/utils/atomsbase_adapter.jl`) は
`AbstractSystem → TOML 形状の dict → Config4System` という経路を取った。
kwarg コンストラクタ `SCEBasis(; lattice, kd, ...)` も同様に
`Config4System` を経由する見込み。

`Config4System` は元々 TOML パース専用の中間型で、validation と
OffsetArray 構築のロジックを持つ。新 API ではすべての入力経路が
これを単一の parse / validation ターゲットとして共有する形になっている。

## 検討事項

理想の実装にとって `Config4System` は必要か、不要か。

- **残す案**: 全入力経路（TOML / dict / kwarg / AtomsBase）の単一の
  parse・validation ターゲット。検証ロジックが 1 箇所に集約され、
  入力経路が増えても検証は共通。
- **バイパス案**: `Config4System` を介さず、各コンストラクタが
  `Structure` + interaction データを直接組み立てる。中間型が減り、
  「TOML 形状の dict」という人工的な中間表現がなくなる。ただし
  validation の集約先を別途設計する必要がある。

## スコープ

別案件。SCE 公開 API リファクタ spec の実装が完了し、4 型 API が
固まった後に着手判断する。
