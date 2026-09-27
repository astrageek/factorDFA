#' Детерминированный факторный анализ
#'
#' Единая точка входа: выбирает метод по аргументу `method`.
#'
#' @inheritParams chain_substitution
#' @param method Метод: `"chain"` — цепные подстановки, `"abs"` — абсолютные
#'   разницы, `"rel"` — относительные разницы, `"index"` — индексный,
#'   `"integral"` — интегральный.
#' @return Объект класса `dfa`.
#' @seealso [dfa_compare()] — все методы сразу.
#' @export
#' @examples
#' dfa(V ~ N * W, base = c(N = 100, W = 10), actual = c(N = 120, W = 12),
#'     method = "integral")
dfa <- function(model, base, actual = NULL,
                method = c("chain", "abs", "rel", "index", "integral"),
                order = NULL) {
  method <- match.arg(method)
  fn <- switch(method,
    chain = chain_substitution,
    abs = abs_diff,
    rel = rel_diff,
    index = index_method,
    integral = integral_method
  )
  fn(model, base, actual, order = order)
}

#' Сравнение методов факторного анализа
#'
#' Считает влияние факторов всеми подходящими методами и сводит результаты
#' в одну таблицу. Методы, которые не подходят для модели, пропускаются
#' с объяснением причины.
#'
#' @inheritParams dfa
#' @param methods Какие методы сравнивать.
#' @return `data.frame` класса `dfa_compare`: строки — факторы и итог,
#'   колонки — методы. Причины пропуска — в атрибуте `"skipped"`.
#' @export
#' @examples
#' dfa_compare(V ~ N * W, base = c(N = 100, W = 10), actual = c(N = 120, W = 12))
dfa_compare <- function(model, base, actual = NULL, order = NULL,
                        methods = c("chain", "abs", "rel", "index", "integral")) {
  m <- .parse_model(model)
  order <- .resolve_order(m, order)
  tab <- data.frame(factor = c(order, "Итого"), stringsAsFactors = FALSE)
  skipped <- character()
  for (mt in methods) {
    res <- tryCatch(dfa(model, base, actual, method = mt, order = order),
                    error = function(e) e)
    if (inherits(res, "error")) {
      skipped[.method_labels[[mt]]] <- conditionMessage(res)
      next
    }
    eff <- res$effects$effect[match(order, res$effects$factor)]
    tab[[.method_labels[[mt]]]] <- c(eff, sum(eff))
  }
  structure(tab, class = c("dfa_compare", "data.frame"), skipped = skipped)
}

#' @export
print.dfa_compare <- function(x, digits = 2, ...) {
  out <- x
  class(out) <- "data.frame"
  names(out)[1L] <- "Фактор"
  for (j in seq_along(out)[-1L]) out[[j]] <- .fmt_signed(out[[j]], digits)
  cat("Сравнение методов факторного анализа\n\n")
  print(out, row.names = FALSE, right = FALSE)
  sk <- attr(x, "skipped")
  if (length(sk)) {
    cat("\nНе применимы к этой модели:\n")
    for (n in names(sk)) cat(" - ", n, ": ", sk[[n]], "\n", sep = "")
  }
  invisible(x)
}
