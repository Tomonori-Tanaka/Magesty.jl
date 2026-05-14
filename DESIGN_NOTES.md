# Design Notes Index

設計検討・調査結果・バックログのインデックス。詳細は各リンク先を参照。
着手済みの開発単位は `docs/specs/` を参照。実装済みの履歴は `.claude/bench_log.md` と `git log`。

ディレクトリ運用ルールは [`docs/design-notes/README.md`](docs/design-notes/README.md)。

## 設計提案

| トピック | Status | 最終更新 |
|---|---|---|
| [SCE 公開 API: 4 型構成 + StatsAPI 化](docs/design-notes/sce-public-api.md) | 未着手 | 2026-05-14 |
| [自作コンテナ (SortedContainer / CountingContainer) を DataStructures.jl で置き換える](docs/design-notes/replace-sorted-container.md) | CountingContainer 完了 / SortedContainer 未着手 | 2026-05-14 |
| [Optimize.jl estimator dispatch](docs/design-notes/estimator-dispatch.md) | 完了 (branch `refactor/estimator-dispatch`) | 2026-05-13 |

## リファクタリング進捗

- [リファクタリング候補スイープ R1-R11](docs/design-notes/refactor-sweep.md) — R1-R7, R9-R11 完了。R8 は Plan B 完了、Plan C は `replace-sorted-container.md` で吸収予定 (2026-05-14)

## 調査結果

- [SpheriCart.jl 採用の可否](docs/design-notes/investigations/sphericart-adoption.md) — 結論: 不採用（2-3× 遅い）(2026-05-11)

## パフォーマンスバックログ

- [Zₗₘ / Legendre / SH バッファ系の改善案](docs/design-notes/backlog.md) — 軽量な思いつきメモ、保留候補
