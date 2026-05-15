# Design Notes Index

設計検討・調査結果・バックログのインデックス。詳細は各リンク先を参照。
着手済みの開発単位は `docs/specs/` を参照。実装済みの履歴は `.claude/bench_log.md` と `git log`。

ディレクトリ運用ルールは [`docs/design-notes/README.md`](docs/design-notes/README.md)。

## 設計提案

| トピック | Status | 最終更新 |
|---|---|---|
| [Config4System の役割の見直し](docs/design-notes/config4system-role.md) | 未着手（SCE 公開 API spec 完了後に判断） | 2026-05-14 |

完了済みの設計提案は対応する spec (`docs/specs/`) に統合されているため本インデックスから外している。spec フォルダ自身が公式の履歴であり、design-note 本体は削除している（commit history で参照可能）。

## リファクタリング進捗

- [リファクタリング候補スイープ R1-R11](docs/design-notes/refactor-sweep.md) — R1-R7, R9-R11 完了。R8 は Plan B 完了、Plan C はソートコンテナ置換リファクタで吸収済み (2026-05-14)

## 調査結果

- [SpheriCart.jl 採用の可否](docs/design-notes/investigations/sphericart-adoption.md) — 結論: 不採用（2-3× 遅い）(2026-05-11)

## パフォーマンスバックログ

- [Zₗₘ / Legendre / SH バッファ系の改善案](docs/design-notes/backlog.md) — 軽量な思いつきメモ、保留候補
- [Post-Step 7 cleanup](docs/design-notes/post-step7-cleanup.md) — SCE 公開 API 破壊コミット完了後の追走項目（hot-path 型安定性、XML I/O 再パース、変数名と中身の乖離スイープなど）

