# True ATT of a (semi-)synthetic panel

Convenience accessor for the ground-truth average treatment effect on
the treated of a simulated panel. It reads the per-cell true effect —
the `tau` column produced by
[`sim_semisynthetic()`](https://takuma1102.github.io/cfcompare/reference/sim_semisynthetic.md),
or `y - y0` for a
[`sim_panel()`](https://takuma1102.github.io/cfcompare/reference/sim_panel.md)
panel — and averages it over the actively treated cells, giving the
number estimators are trying to recover. Pair it with
[`panel_compare()`](https://takuma1102.github.io/cfcompare/reference/panel_compare.md)
or
[`panel_rmse()`](https://takuma1102.github.io/cfcompare/reference/panel_rmse.md)
to score estimates against the truth without re-deriving it by hand.

## Usage

``` r
true_att(data, treatment = "w", tau = "tau")
```

## Arguments

- data:

  A long `data.frame` with a treatment indicator and either a `tau`
  column (true per-cell effect) or both `y` and `y0` columns.

- treatment:

  Name of the 0/1 treatment column (default `"w"`).

- tau:

  Name of the true-effect column (default `"tau"`). When that column is
  absent the per-cell effect is taken as `y - y0`.

## Value

A single numeric: the true ATT averaged over treated cells.

## See also

[`sim_semisynthetic()`](https://takuma1102.github.io/cfcompare/reference/sim_semisynthetic.md),
[`sim_panel()`](https://takuma1102.github.io/cfcompare/reference/sim_panel.md),
[`panel_compare()`](https://takuma1102.github.io/cfcompare/reference/panel_compare.md)

## Examples

``` r
real <- sim_panel(N = 16, T = 10, n_treated = 0L, att = 0, seed = 1)
ss <- sim_semisynthetic(real, "y", "id", "t",
                        n_treated = 6, t0 = 8, att = 3, seed = 2)
true_att(ss)   # 3
#> [1] 3
```
