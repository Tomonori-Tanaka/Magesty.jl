# スタイルガイド

Julia 公式スタイルガイドおよび DFTK スタイルガイドに準拠する。
基本方針: **一貫性よりも読みやすさを優先**。

## 命名規則

| 対象 | 規則 | 例 |
|------|------|----|
| モジュール・型 | `UpperCamelCase` | `SpinCluster`, `BasisSet` |
| 関数・変数 | `lowercase` / `snake_case` | `calc_energy`, `num_atoms` |
| 定数 | `UPPER_SNAKE_CASE` | `MAX_ITER` |
| 内部ヘルパー関数 | 先頭に `_` | `_compute_aux` |
| 破壊的関数 | 末尾に `!` | `update_spin!` |

## `for` ループ

- **インデックス（数値カウンタ）** には `=` を使う
- **コレクションの要素** には `in` を使う

```julia
# インデックスループ → =
for i = 1:num_atoms
for isym = 1:n_operations

# 要素ループ → in
for atom in atoms
for (salc_idx, key_group) in enumerate(salc_list)
```

この区別により、変数がインデックスか値かが読んだだけで分かる。

## 関数引数の順序

Julia 公式の推奨順序に従う:

```
function → IO → 変更される引数 → 型 → 変更されない引数 → キーワード引数
```

## 型アノテーション

- 公開 API の引数には型アノテーションを付ける
- 戻り値の型も明示する（`::Float64` 等）
- `where` には必ず波括弧を付ける

```julia
# Good
function calc_energy(sc::SpinCluster, spin_config::AbstractMatrix{<:Real})::Float64
function foo(x::T) where {T <: AbstractFloat}

# Bad
function foo(x::T) where T <: AbstractFloat
```

## 名前付きタプル

明示的な構文 `(; var=val)` を使う:

```julia
# Good
return (; energy=E, forces=F)
(; energy, forces) = result

# Bad
return (energy=E, forces=F)
```

## その他

- インデント: スペース 4 つ
- 行長: 最大 92 文字
- `if`/`while` の条件式を括弧でくくらない
- 型チェックには `isa` / `<:` を使う（`==` は使わない）
- 空のコールバックには `identity` を渡す
- 型パイラシー（他のパッケージが定義した型への Base メソッド追加）を避ける
