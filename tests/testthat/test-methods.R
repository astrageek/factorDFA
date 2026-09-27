b2 <- c(a = 100, b = 10)
a2 <- c(a = 120, b = 12)

eff <- function(res) stats::setNames(res$effects$effect, res$effects$factor)

test_that("цепные подстановки: y = a * b", {
  r <- chain_substitution(y ~ a * b, b2, a2)
  expect_equal(eff(r), c(a = 200, b = 240))
  expect_equal(r$delta, 440)
  expect_equal(r$residual, 0)
  expect_equal(r$details$steps$y, c(1000, 1200, 1440))
})

test_that("порядок подстановки меняет результат цепного метода", {
  r <- chain_substitution(y ~ a * b, b2, a2, order = c("b", "a"))
  expect_equal(eff(r), c(b = 200, a = 240))
})

test_that("абсолютные разницы совпадают с цепными подстановками", {
  expect_equal(eff(abs_diff(y ~ a * b, b2, a2)), c(a = 200, b = 240))
  r <- abs_diff("P = V * (price - cost)",
                c(V = 1000, price = 50, cost = 35),
                c(V = 1100, price = 52, cost = 38))
  expect_equal(eff(r), c(V = 1500, price = 2200, cost = -3300))
  expect_equal(r$delta, 400)
  expect_equal(r$effects$formula[1], "ΔV × (price0 - cost0)")
})

test_that("аддитивная модель", {
  r <- abs_diff(y ~ a + b - c, c(a = 10, b = 5, c = 3), c(a = 12, b = 4, c = 6))
  expect_equal(eff(r), c(a = 2, b = -1, c = -3))
})

test_that("относительные разницы", {
  r <- rel_diff(y ~ a * b, b2, a2)
  expect_equal(eff(r), c(a = 200, b = 240))
  x0 <- c(V = 2000, UR = 3, CM = 40); x1 <- c(V = 2200, UR = 2.8, CM = 45)
  expect_equal(eff(rel_diff(MZ ~ V * UR * CM, x0, x1)),
               eff(chain_substitution(MZ ~ V * UR * CM, x0, x1)))
  expect_error(rel_diff(y ~ a / b, b2, a2), "относительные разницы")
})

test_that("индексный метод, в том числе агрегатные индексы", {
  r <- index_method(y ~ a * b, b2, a2)
  expect_equal(r$details$indices$index, c(1.2, 1.2))
  expect_equal(r$details$total_index, 1.44)
  r2 <- index_method(R ~ q * p,
                     data.frame(q = c(100, 50), p = c(20, 40)),
                     data.frame(q = c(110, 45), p = c(22, 44)))
  expect_equal(c(r2$y0, r2$y1), c(4000, 4400))
  expect_equal(eff(r2), c(q = 0, p = 400))
  expect_equal(r2$details$indices$index, c(1, 1.1))
  expect_error(index_method(y ~ a + b, b2, a2), "индексный")
})

test_that("интегральный метод: произведение и частное", {
  expect_equal(eff(integral_method(y ~ a * b, b2, a2)), c(a = 220, b = 220))
  r <- integral_method(W ~ V / N, c(V = 1000, N = 40), c(V = 1500, N = 50))
  ev <- 500 / 10 * log(50 / 40)
  expect_equal(eff(r), c(V = ev, N = 5 - ev))
  expect_lt(abs(r$residual), 1e-9)
})

test_that("интегральный метод: три фактора по учебной формуле", {
  x0 <- c(a = 10, b = 4, c = 2); x1 <- c(a = 12, b = 5, c = 3)
  d <- x1 - x0
  r <- integral_method(y ~ a * b * c, x0, x1)
  ea <- 0.5 * d[["a"]] * (x0[["b"]] * x1[["c"]] + x1[["b"]] * x0[["c"]]) + prod(d) / 3
  expect_equal(eff(r)[["a"]], ea)
  expect_lt(abs(r$residual), 1e-9)
  r2 <- integral_method(y ~ a * b * c, x0, x1, order = c("c", "b", "a"))
  expect_equal(eff(r2), eff(r)[c("c", "b", "a")])
})

test_that("dfa() и dfa_compare()", {
  expect_equal(eff(dfa(y ~ a * b, b2, a2, method = "integral")), c(a = 220, b = 220))
  cmp <- dfa_compare(W ~ V / N, c(V = 1000, N = 40), c(V = 1500, N = 50))
  expect_true(all(c("цепные подстановки", "индексный", "интегральный") %in% names(cmp)))
  expect_true("абсолютные разницы" %in% names(attr(cmp, "skipped")))
  expect_output(print(cmp), "Сравнение методов")
})

test_that("печать и график работают", {
  for (mt in c("chain", "abs", "rel", "index", "integral")) {
    r <- dfa(V ~ N * W, c(N = 100, W = 10), c(N = 120, W = 12), method = mt)
    expect_output(print(r), "Баланс сходится")
  }
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_silent(plot(r))
})
