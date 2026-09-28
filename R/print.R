#' Печать результата факторного анализа
#'
#' @param x Объект класса `dfa`.
#' @param digits Знаков после запятой.
#' @param details Показывать ли расчётные таблицы метода.
#' @param ... Не используется.
#' @export
print.dfa <- function(x, digits = 2, details = TRUE, ...) {
  lhs <- x$model$lhs
  cat("Детерминированный факторный анализ\n")
  cat("Метод: ", .method_labels[[x$method]], "\n", sep = "")
  cat("Модель: ", lhs, " = ", x$model$text, "\n", sep = "")
  if (x$n_items > 1L) {
    cat("Объектов: ", x$n_items, " (показатель суммируется)\n", sep = "")
  }
  if (!x$method %in% c("integral", "log", "shapley")) {
    cat("Порядок подстановки: ", paste(x$order, collapse = " → "), "\n", sep = "")
  }
  cat(lhs, "0 = ", .fmt(x$y0, digits), ";  ", lhs, "1 = ", .fmt(x$y1, digits), "\n", sep = "")
  pct <- if (x$y0 != 0) paste0(" (", .fmt_signed(x$delta / abs(x$y0) * 100, digits), "%)") else ""
  cat("Δ", lhs, " = ", .fmt_signed(x$delta, digits), pct, "\n\n", sep = "")

  eff <- x$effects
  tab <- data.frame(
    `Фактор` = c(eff$factor, "Итого"),
    `Влияние` = .fmt_signed(c(eff$effect, sum(eff$effect)), digits),
    `Доля, %` = .fmt(c(eff$share, if (x$delta != 0) sum(eff$share) else NA), digits),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  if (!is.null(eff$formula)) tab[["Формула"]] <- c(eff$formula, "")
  cat("Влияние факторов:\n")
  print(tab, row.names = FALSE, right = FALSE)

  tol <- 1e-8 * max(1, abs(x$delta))
  if (abs(x$residual) <= tol) {
    cat("\nБаланс сходится: сумма влияний = Δ", lhs, "\n", sep = "")
  } else {
    cat("\nВнимание: баланс не сходится, расхождение ", .fmt(x$residual, 6), "\n", sep = "")
  }

  if (details) .print_details(x, digits)
  invisible(x)
}

.print_details <- function(x, digits) {
  d <- x$details
  num_cols <- function(df, cols, signed = character()) {
    for (cl in cols) df[[cl]] <- if (cl %in% signed) .fmt_signed(df[[cl]], digits) else .fmt(df[[cl]], digits)
    df
  }
  if (x$method == "chain") {
    tbl <- d$steps
    numeric_cols <- names(tbl)[vapply(tbl, is.numeric, logical(1))]
    tbl <- num_cols(tbl, numeric_cols, signed = "Влияние")
    cat("\nРасчётная таблица подстановок:\n")
    print(tbl, row.names = FALSE, right = FALSE)
  } else if (x$method == "abs") {
    tbl <- num_cols(d$deltas, c("base", "actual", "delta"), signed = "delta")
    names(tbl) <- c("Фактор", "База", "Отчёт", "Δ")
    cat("\nАбсолютные отклонения факторов:\n")
    print(tbl, row.names = FALSE, right = FALSE)
  } else if (x$method == "rel") {
    tbl <- num_cols(d$changes, c("base", "actual", "change_pct"), signed = "change_pct")
    names(tbl) <- c("Фактор", "База", "Отчёт", "Δ%")
    cat("\nОтносительные отклонения факторов:\n")
    print(tbl, row.names = FALSE, right = FALSE)
  } else if (x$method == "index") {
    tbl <- d$indices
    tbl$index <- formatC(tbl$index, format = "f", digits = 4)
    tbl <- num_cols(tbl, c("change_pct", "numerator", "denominator"), signed = "change_pct")
    names(tbl) <- c("Фактор", "Индекс", "Изм., %", "Числитель", "Знаменатель", "Формула")
    cat("\nФакторные индексы:\n")
    print(tbl, row.names = FALSE, right = FALSE)
    cat("Общий индекс I(", x$model$lhs, ") = ", formatC(d$total_index, format = "f", digits = 4),
        " = произведение факторных индексов\n", sep = "")
  } else if (x$method == "log") {
    tbl <- d$indices
    tbl$index <- formatC(tbl$index, format = "f", digits = 4)
    tbl$log_index <- formatC(tbl$log_index, format = "f", digits = 6)
    names(tbl) <- c("Фактор", "Индекс", "ln индекса")
    cat("\nФакторные индексы (каждый фактор меняется отдельно от базы):\n")
    print(tbl, row.names = FALSE, right = FALSE)
    cat("Общий индекс I(", x$model$lhs, ") = ", formatC(d$total_index, format = "f", digits = 4),
        "; изменение распределено пропорционально логарифмам индексов\n", sep = "")
  } else if (x$method == "shapley") {
    cat("\nВлияние усреднено по всем ", d$n_orders,
        " порядкам подстановки; результат не зависит от порядка факторов.\n", sep = "")
  } else if (x$method == "integral") {
    cat("\nНеразложимый остаток распределён между факторами; результат не зависит",
        "от порядка факторов.\n")
  }
}

#' Таблица влияния факторов
#'
#' @param x Объект класса `dfa`.
#' @param ... Не используется.
#' @return `data.frame` с колонками `factor`, `effect`, `share`, `formula`.
#' @export
as.data.frame.dfa <- function(x, ...) x$effects

#' Каскадная диаграмма (водопад) влияния факторов
#'
#' @param x Объект класса `dfa`.
#' @param col_up,col_down,col_total Цвета роста, снижения и итоговых столбцов.
#' @param main Заголовок.
#' @param digits Знаков после запятой в подписях.
#' @param ... Передаётся в [graphics::plot.window()].
#' @export
#' @examples
#' res <- chain_substitution(V ~ N * W, c(N = 100, W = 10), c(N = 120, W = 12))
#' plot(res)
plot.dfa <- function(x, col_up = "#2E7D32", col_down = "#C62828", col_total = "#546E7A",
                     main = NULL, digits = 2, ...) {
  lhs <- x$model$lhs
  eff <- x$effects$effect
  labs <- c(paste0(lhs, "0"), x$effects$factor, paste0(lhs, "1"))
  n <- length(labs)
  ends <- c(x$y0, x$y0 + cumsum(eff), x$y1)
  starts <- c(0, x$y0 + c(0, cumsum(eff))[seq_along(eff)], 0)
  cols <- c(col_total, ifelse(eff >= 0, col_up, col_down), col_total)
  lo <- pmin(starts, ends)
  hi <- pmax(starts, ends)
  yr <- range(c(lo, hi, 0))
  yr[2L] <- yr[2L] + diff(yr) * 0.1
  graphics::plot.new()
  graphics::plot.window(xlim = c(0.4, n + 0.6), ylim = yr, ...)
  graphics::rect(seq_len(n) - 0.35, lo, seq_len(n) + 0.35, hi, col = cols, border = NA)
  graphics::segments(seq_len(n - 1L) + 0.35, ends[-n], seq_len(n - 1L) + 0.65, ends[-n], lty = 3)
  graphics::abline(h = 0, col = "grey60")
  graphics::axis(1, at = seq_len(n), labels = labs, tick = FALSE)
  graphics::axis(2, las = 1)
  graphics::text(seq_len(n), hi, c(.fmt(x$y0, digits), .fmt_signed(eff, digits), .fmt(x$y1, digits)),
                 pos = 3, cex = 0.8, xpd = TRUE)
  if (is.null(main)) main <- paste0("Факторный анализ ", lhs, ": ", .method_labels[[x$method]])
  graphics::title(main = main)
  invisible(x)
}
