#' Метод цепных подстановок
#'
#' Последовательно заменяет базисные значения факторов фактическими и
#' считает влияние каждого фактора как разность соседних условных
#' значений результативного показателя. Подходит для любых моделей:
#' мультипликативных, аддитивных, кратных и смешанных.
#'
#' @param model Модель: формула (`V ~ N * W`) или строка (`"V = N * W"`).
#'   Если левая часть не указана, показатель называется `y`.
#' @param base Базисные (плановые) значения факторов: именованный вектор,
#'   список или `data.frame` (строки — объекты, колонки — факторы).
#'   Можно передать одну таблицу без `actual`: с колонками
#'   `factor` / `base` / `actual` (или `фактор` / `план` / `факт`),
#'   либо из двух строк (база и отчёт), где факторы — колонки.
#' @param actual Фактические (отчётные) значения в том же формате, что `base`.
#' @param order Порядок подстановки факторов. По умолчанию — порядок
#'   появления в модели. Обычно сначала количественные факторы, потом
#'   качественные.
#'
#' @return Объект класса `dfa`: влияние факторов (`$effects`),
#'   значения `$y0`, `$y1`, `$delta` и расчётная таблица (`$details$steps`).
#' @export
#' @examples
#' # Выручка = численность * выработка
#' chain_substitution(V ~ N * W,
#'                    base   = c(N = 100, W = 10),
#'                    actual = c(N = 120, W = 12))
chain_substitution <- function(model, base, actual = NULL, order = NULL) {
  m <- .parse_model(model)
  order <- .resolve_order(m, order)
  d <- .parse_data(base, actual, m$vars)
  steps <- .chain_steps(m, d, order)
  effects <- stats::setNames(diff(steps$y), order)
  formulas <- vapply(seq_along(order), function(i) {
    paste0(m$lhs, if (i == length(order)) "1" else paste0("_усл", i), " − ",
           m$lhs, if (i == 1L) "0" else paste0("_усл", i - 1L))
  }, character(1))
  .new_dfa("chain", m, d, order, effects,
           details = list(steps = steps$table),
           formulas = stats::setNames(formulas, order))
}

#' Метод абсолютных разниц
#'
#' Влияние фактора равно его абсолютному приросту, умноженному на
#' значения остальных факторов: предшествующих — фактических, последующих —
#' базисных. Применяется к мультипликативным (`a * b * c`), аддитивным
#' (`a + b - c`) и смешанным (`a * (b - c)`) моделям, где каждый фактор
#' встречается один раз и нет деления.
#'
#' @inheritParams chain_substitution
#' @return Объект класса `dfa`.
#' @export
#' @examples
#' # Прибыль = объём * (цена - себестоимость)
#' abs_diff("P = V * (price - cost)",
#'          base   = c(V = 1000, price = 50, cost = 35),
#'          actual = c(V = 1100, price = 52, cost = 38))
abs_diff <- function(model, base, actual = NULL, order = NULL) {
  m <- .parse_model(model)
  .check_structure(m, c("*", "+", "-", "("), "абсолютные разницы",
                   "Используйте chain_substitution() или integral_method().")
  order <- .resolve_order(m, order)
  d <- .parse_data(base, actual, m$vars)
  effects <- stats::setNames(numeric(length(order)), order)
  formulas <- stats::setNames(character(length(order)), order)
  for (i in seq_along(order)) {
    v <- order[i]
    before <- order[seq_len(i - 1L)]
    after <- order[-seq_len(i)]
    point <- d$base
    for (u in before) point[[u]] <- d$actual[[u]]
    dv <- stats::D(m$expr, v)
    while (is.call(dv) && identical(dv[[1L]], as.name("("))) dv <- dv[[2L]]
    effects[v] <- sum((d$actual[[v]] - d$base[[v]]) * .eval_expr(dv, point))
    slope <- .sub_text(dv, before, after)
    formulas[v] <- switch(slope,
      "1" = paste0("Δ", v),
      "-1" = paste0("−Δ", v),
      paste0("Δ", v, " × (", slope, ")")
    )
  }
  deltas <- data.frame(
    factor = order,
    base = vapply(order, function(v) sum(d$base[[v]]), numeric(1)),
    actual = vapply(order, function(v) sum(d$actual[[v]]), numeric(1)),
    stringsAsFactors = FALSE, row.names = NULL
  )
  deltas$delta <- deltas$actual - deltas$base
  .new_dfa("abs", m, d, order, effects,
           details = list(deltas = deltas), formulas = formulas)
}

#' Метод относительных разниц
#'
#' Влияние первого фактора равно базисному значению показателя, умноженному
#' на относительный прирост фактора (в процентах) / 100; влияние каждого
#' следующего — показателю с учётом влияния предыдущих факторов,
#' умноженному на прирост этого фактора. Применяется только к чисто
#' мультипликативным моделям (`a * b * c`) для одного объекта.
#'
#' @inheritParams chain_substitution
#' @return Объект класса `dfa`.
#' @export
#' @examples
#' # Материальные затраты = объём * норма расхода * цена материала
#' rel_diff(MZ ~ V * UR * CM,
#'          base   = c(V = 2000, UR = 3, CM = 40),
#'          actual = c(V = 2200, UR = 2.8, CM = 45))
rel_diff <- function(model, base, actual = NULL, order = NULL) {
  m <- .parse_model(model)
  .check_structure(m, c("*", "("), "относительные разницы",
                   "Метод работает только для произведения факторов (a * b * c).")
  order <- .resolve_order(m, order)
  d <- .parse_data(base, actual, m$vars)
  if (nrow(d$base) > 1L) {
    stop("Метод относительных разниц применяется к одному объекту. ",
         "Для нескольких объектов используйте chain_substitution() или index_method().",
         call. = FALSE)
  }
  zero <- order[vapply(order, function(v) d$base[[v]] == 0, logical(1))]
  if (length(zero)) {
    stop("Базисное значение равно нулю, относительный прирост не определён: ",
         paste(zero, collapse = ", "), call. = FALSE)
  }
  pct <- vapply(order, function(v) (d$actual[[v]] - d$base[[v]]) / d$base[[v]] * 100, numeric(1))
  cum <- .total(m, d$base)
  effects <- stats::setNames(numeric(length(order)), order)
  formulas <- stats::setNames(character(length(order)), order)
  for (i in seq_along(order)) {
    v <- order[i]
    effects[v] <- cum * pct[[v]] / 100
    cum <- cum + effects[v]
    prev <- order[seq_len(i - 1L)]
    formulas[v] <- if (!length(prev)) {
      paste0(m$lhs, "0 × Δ%", v, " / 100")
    } else {
      paste0("(", m$lhs, "0 + ", paste0("Δ", m$lhs, "(", prev, ")", collapse = " + "),
             ") × Δ%", v, " / 100")
    }
  }
  pct_tbl <- data.frame(
    factor = order,
    base = vapply(order, function(v) d$base[[v]], numeric(1)),
    actual = vapply(order, function(v) d$actual[[v]], numeric(1)),
    change_pct = unname(pct),
    stringsAsFactors = FALSE, row.names = NULL
  )
  .new_dfa("rel", m, d, order, effects,
           details = list(changes = pct_tbl), formulas = formulas)
}

#' Индексный метод
#'
#' Строит факторные индексы как отношения соседних условных значений
#' показателя (агрегатные индексы при нескольких объектах, например
#' \eqn{I_q = \sum q_1 p_0 / \sum q_0 p_0}). Произведение факторных индексов
#' равно общему индексу, а разности числителя и знаменателя дают абсолютное
#' влияние факторов. Применяется к мультипликативным и кратным моделям
#' (`a * b`, `a / b`, `a * b / c`).
#'
#' @inheritParams chain_substitution
#' @return Объект класса `dfa`; индексы — в `$details$indices`,
#'   общий индекс — в `$details$total_index`.
#' @export
#' @examples
#' # Выручка по двум товарам: количество * цена
#' index_method(R ~ q * p,
#'              base   = data.frame(q = c(100, 50), p = c(20, 40)),
#'              actual = data.frame(q = c(110, 45), p = c(22, 44)))
index_method <- function(model, base, actual = NULL, order = NULL) {
  m <- .parse_model(model)
  .check_structure(m, c("*", "/", "("), "индексный",
                   "Метод работает для мультипликативных и кратных моделей (a * b, a / b).")
  order <- .resolve_order(m, order)
  d <- .parse_data(base, actual, m$vars)
  y <- .chain_steps(m, d, order)$y
  if (any(y[-length(y)] == 0)) {
    stop("Одно из условных значений показателя равно нулю, индекс не определён.", call. = FALSE)
  }
  idx <- y[-1L] / y[-length(y)]
  sigma <- if (nrow(d$base) > 1L) "Σ" else ""
  texts <- vapply(0:length(order), function(i) {
    paste0(sigma, .sub_text(m$expr, order[seq_len(i)], order[-seq_len(i)]))
  }, character(1))
  num <- texts[-1L]
  den <- texts[-length(texts)]
  indices <- data.frame(
    factor = order, index = idx, change_pct = (idx - 1) * 100,
    numerator = y[-1L], denominator = y[-length(y)],
    formula = paste0("(", num, ") / (", den, ")"),
    stringsAsFactors = FALSE
  )
  .new_dfa("index", m, d, order, stats::setNames(diff(y), order),
           details = list(indices = indices, total_index = y[length(y)] / y[1L]),
           formulas = stats::setNames(paste0(num, " − ", den), order))
}

#' Интегральный метод
#'
#' Распределяет неразложимый остаток между факторами без зависимости от
#' порядка подстановки. Влияние фактора \eqn{x_i} равно интегралу
#' \eqn{\int_0^1 \partial f / \partial x_i \cdot \Delta x_i \, dt} вдоль
#' прямой от базисных значений к фактическим; интеграл вычисляется
#' квадратурой Гаусса–Лежандра (для многочленов — точно). Работает для
#' мультипликативных, кратных и смешанных моделей. Для типовых моделей
#' выводятся учебные формулы, например для \eqn{y = ab}:
#' \eqn{\Delta y_a = \Delta a \, b_0 + \frac12 \Delta a \Delta b}.
#'
#' @inheritParams chain_substitution
#' @param order Порядок вывода факторов в таблице (на результат не влияет).
#' @param n_nodes Число узлов квадратуры.
#' @return Объект класса `dfa`.
#' @export
#' @examples
#' # Средняя выработка = выручка / численность
#' integral_method(W ~ V / N,
#'                 base   = c(V = 1000, N = 40),
#'                 actual = c(V = 1500, N = 50))
integral_method <- function(model, base, actual = NULL, order = NULL, n_nodes = 40L) {
  m <- .parse_model(model)
  order <- .resolve_order(m, order)
  d <- .parse_data(base, actual, m$vars)
  gl <- .gauss_legendre(n_nodes)
  points <- lapply(gl$t, function(t) {
    p <- d$base
    for (u in m$vars) p[[u]] <- d$base[[u]] + t * (d$actual[[u]] - d$base[[u]])
    p
  })
  effects <- stats::setNames(numeric(length(m$vars)), m$vars)
  for (v in m$vars) {
    dx <- d$actual[[v]] - d$base[[v]]
    if (all(dx == 0)) next
    dv <- tryCatch(stats::D(m$expr, v), error = function(e) {
      stop(sprintf("Не удалось продифференцировать модель по %s: %s", v, conditionMessage(e)),
           call. = FALSE)
    })
    vals <- vapply(points, function(p) sum(.eval_expr(dv, p) * dx), numeric(1))
    if (any(!is.finite(vals))) {
      stop(sprintf("Производная по %s не определена между базой и отчётом ", v),
           "(например, знаменатель проходит через ноль).", call. = FALSE)
    }
    effects[v] <- sum(gl$w * vals)
  }
  .new_dfa("integral", m, d, order, effects,
           details = list(n_nodes = n_nodes),
           formulas = .integral_formulas(m, order))
}

.flatten_product <- function(e) {
  if (is.name(e)) return(list(e))
  if (is.call(e) && identical(e[[1L]], as.name("("))) return(.flatten_product(e[[2L]]))
  if (is.call(e) && identical(e[[1L]], as.name("*")) && length(e) == 3L) {
    l <- .flatten_product(e[[2L]])
    r <- .flatten_product(e[[3L]])
    if (is.null(l) || is.null(r)) return(NULL)
    return(c(l, r))
  }
  NULL
}

.integral_formulas <- function(m, order) {
  e <- m$expr
  prod_vars <- vapply(.flatten_product(e), as.character, character(1))
  out <- stats::setNames(paste0("∫ ∂", m$lhs, "/∂", order, " · Δ", order, " dt"), order)
  if (length(prod_vars) == 2L && !anyDuplicated(prod_vars)) {
    a <- prod_vars[1L]; b <- prod_vars[2L]
    out[a] <- paste0("Δ", a, " × ", b, "0 + ½ × Δ", a, " × Δ", b)
    out[b] <- paste0("Δ", b, " × ", a, "0 + ½ × Δ", a, " × Δ", b)
  } else if (length(prod_vars) == 3L && !anyDuplicated(prod_vars)) {
    for (i in 1:3) {
      x <- prod_vars[i]
      o <- prod_vars[-i]
      out[x] <- paste0("½ × Δ", x, " × (", o[1], "0 × ", o[2], "1 + ", o[1], "1 × ", o[2],
                       "0) + ⅓ × Δ", prod_vars[1], " × Δ", prod_vars[2], " × Δ", prod_vars[3])
    }
  } else if (is.call(e) && identical(e[[1L]], as.name("/")) && length(e) == 3L &&
             is.name(e[[2L]]) && is.name(e[[3L]])) {
    a <- as.character(e[[2L]]); b <- as.character(e[[3L]])
    out[a] <- paste0("Δ", a, " / Δ", b, " × ln(", b, "1 / ", b, "0)")
    out[b] <- paste0("Δ", m$lhs, " − Δ", m$lhs, "(", a, ")")
  }
  out
}
