eff <- function(res) stats::setNames(res$effects$effect, res$effects$factor)
expected <- c(N = 200, W = 240)

test_that("таблица factor / base / actual", {
  df <- data.frame(factor = c("N", "W"), base = c(100, 10), actual = c(120, 12))
  expect_equal(eff(chain_substitution(V ~ N * W, df)), expected)
})

test_that("таблица с русскими названиями колонок", {
  df <- data.frame(`Фактор` = c("N", "W"), `План` = c(100, 10), `Факт` = c(120, 12),
                   check.names = FALSE)
  expect_equal(eff(chain_substitution(V ~ N * W, df)), expected)
})

test_that("таблица из двух строк: база и отчёт", {
  df <- data.frame(period = c("2024", "2025"), N = c(100, 120), W = c(10, 12))
  expect_equal(eff(chain_substitution(V ~ N * W, df)), expected)
})

test_that("списки и строковая модель", {
  r <- chain_substitution("V = N * W", list(N = 100, W = 10), list(N = 120, W = 12))
  expect_equal(eff(r), expected)
  expect_equal(r$model$lhs, "V")
})

test_that("понятные ошибки", {
  expect_error(chain_substitution(V ~ N * W, c(N = 100), c(N = 120, W = 12)), "нет значений")
  expect_error(chain_substitution(V ~ N * W, c(N = 100, W = 10), c(N = 120, W = 12),
                                  order = c("N")), "order")
  expect_error(chain_substitution(V ~ N * W, c(N = 100, W = NA), c(N = 120, W = 12)), "NA")
})
