# Design Notes Index

設計検討・調査結果・バックログのインデックス。詳細は各リンク先を参照。
着手済みの開発単位は `docs/specs/` を参照。実装済みの履歴は `.claude/bench_log.md` と `git log`。

ディレクトリ運用ルールは [`docs/design-notes/README.md`](docs/design-notes/README.md)。

## 設計提案

| トピック | Status | 最終更新 |
|---|---|---|
| [LASSO / Adaptive LASSO / Adaptive Ridge estimator の導入](docs/design-notes/lasso-adaptive-estimators.md) | 未着手（GLMNet.jl 採用方針で合意、spec 化前ドラフト） | 2026-05-16 |

完了済みの設計提案は対応する spec (`docs/specs/`) に統合されているため本インデックスから外している。spec フォルダ自身が公式の履歴であり、design-note 本体は削除している（commit history で参照可能）。

## 調査結果

- [SpheriCart.jl 採用の可否](docs/design-notes/investigations/sphericart-adoption.md) — 結論: 不採用（2-3× 遅い）(2026-05-11)

## パフォーマンスバックログ

- [Zₗₘ / Legendre / SH バッファ系の改善案](docs/design-notes/backlog.md) — 軽量な思いつきメモ、保留候補

