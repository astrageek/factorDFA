.method_labels <- c(
  chain    = "цепные подстановки",
  abs      = "абсолютные разницы",
  rel      = "относительные разницы",
  index    = "индексный",
  integral = "интегральный"
)

.parse_model <- function(model) {
  if (inherits(model, "formula")) {
    rhs <- model[[length(model)]]
    lhs <- if (length(model) == 3L) deparse1(model[[2L]]) else "y"
  } else if (is.character(model) && length(model) == 1L) {
    parts <- strsplit(model, "=", fixed = TRUE)[[1L]]
    if (length(parts) == 2L) {
      lhs <- trimws(parts[1L])
      rhs <- str2lang(parts[2L])
    } else if (length(parts) == 1L) {
      lhs <- "y"
      rhs <- str2lang(parts[1L])
    } else {
      stop("Модель должна содержать не больше одного знака '=', например \"V = N * W\".", call. = FALSE)
    }
  } else if (is.call(model) || is.name(model)) {
    lhs <- "y"
    rhs <- model
  } else {
    stop("`model` задаётся формулой (V ~ N * W) или строкой (\"V = N * W\").", call. = FALSE)
  }
  vars <- all.vars(rhs)
  if (!length(vars)) stop("В модели нет ни одного фактора.", call. = FALSE)
  structure(
    list(expr = rhs, lhs = lhs, vars = vars, text = deparse1(rhs)),
    class = "dfa_model"
  )
}

.check_structure <- function(m, allowed, method, hint) {
  nms <- all.names(m$expr)
  bad <- setdiff(setdiff(nms, m$vars), allowed)
  repeated <- m$vars[vapply(m$vars, function(v) sum(nms == v) > 1L, logical(1))]
  if (length(bad) || length(repeated)) {
    stop(sprintf(
      "Метод «%s» не подходит для модели %s = %s: %s. %s",
      method, m$lhs, m$text,
      if (length(bad)) paste0("недопустимые операции ", paste(bad, collapse = " "))
      else paste0("фактор встречается несколько раз: ", paste(repeated, collapse = ", ")),
      hint
    ), call. = FALSE)
  }
}

.resolve_order <- function(m, order) {
  if (is.null(order)) return(m$vars)
  if (!setequal(order, m$vars) || anyDuplicated(order)) {
    stop("`order` должен содержать все факторы модели ровно по одному разу: ",
         paste(m$vars, collapse = ", "), call. = FALSE)
  }
  as.character(order)
}

.aliases <- list(
  factor = c("factor", "фактор", "показатель", "name"),
  base   = c("base", "plan", "план", "базис", "база", "x0", "прошлый"),
  actual = c("actual", "fact", "факт", "отчет", "отчёт", "x1", "текущий")
)

.find_col <- function(nms, key) {
  hit <- which(tolower(trimws(nms)) %in% .aliases[[key]])
  if (length(hit)) hit[1L] else NA_integer_
}

.as_frame <- function(x, what) {
  if (is.data.frame(x)) return(x)
  if (is.list(x) || (is.numeric(x) && !is.null(names(x)))) {
    return(as.data.frame(as.list(x), check.names = FALSE))
  }
  stop(sprintf("`%s` должен быть именованным вектором, списком или data.frame.", what), call. = FALSE)
}

.parse_data <- function(base, actual, vars) {
  if (is.null(actual)) {
    if (!is.data.frame(base)) {
      stop("Передайте `actual` или таблицу с колонками factor / base / actual.", call. = FALSE)
    }
    fc <- .find_col(names(base), "factor")
    bc <- .find_col(names(base), "base")
    ac <- .find_col(names(base), "actual")
    if (!anyNA(c(fc, bc, ac))) {
      f <- as.character(base[[fc]])
      b <- stats::setNames(as.list(as.numeric(base[[bc]])), f)
      a <- stats::setNames(as.list(as.numeric(base[[ac]])), f)
      base <- as.data.frame(b, check.names = FALSE)
      actual <- as.data.frame(a, check.names = FALSE)
    } else if (nrow(base) == 2L) {
      actual <- base[2L, , drop = FALSE]
      base <- base[1L, , drop = FALSE]
    } else {
      stop("Не понял формат таблицы. Нужны колонки factor / base / actual ",
           "или две строки (база и отчёт) с факторами в колонках.", call. = FALSE)
    }
  }
  base <- .as_frame(base, "base")
  actual <- .as_frame(actual, "actual")
  for (p in list(list(base, "base"), list(actual, "actual"))) {
    miss <- setdiff(vars, names(p[[1L]]))
    if (length(miss)) {
      stop(sprintf("В `%s` нет значений для факторов: %s", p[[2L]], paste(miss, collapse = ", ")),
           call. = FALSE)
    }
  }
  if (nrow(base) != nrow(actual)) {
    stop("В `base` и `actual` разное число объектов (строк).", call. = FALSE)
  }
  base <- base[vars]
  actual <- actual[vars]
  for (v in vars) {
    if (!is.numeric(base[[v]]) || !is.numeric(actual[[v]])) {
      stop(sprintf("Фактор %s должен быть числом.", v), call. = FALSE)
    }
    if (anyNA(base[[v]]) || anyNA(actual[[v]])) {
      stop(sprintf("У фактора %s есть пропуски (NA).", v), call. = FALSE)
    }
  }
  rownames(base) <- rownames(actual) <- NULL
  list(base = base, actual = actual)
}

.eval_expr <- function(expr, values) {
  out <- eval(expr, values, baseenv())
  rep_len(as.numeric(out), nrow(values))
}

.total <- function(m, values) sum(.eval_expr(m$expr, values))

.sub_text <- function(expr, actual_vars, base_vars) {
  map <- c(
    stats::setNames(lapply(paste0(actual_vars, "1"), as.name), actual_vars),
    stats::setNames(lapply(paste0(base_vars, "0"), as.name), base_vars)
  )
  gsub(" * ", " × ", deparse1(do.call(substitute, list(expr, map))), fixed = TRUE)
}

.new_dfa <- function(method, m, d, order, effects, details = list(), formulas = NULL) {
  y0 <- .total(m, d$base)
  y1 <- .total(m, d$actual)
  delta <- y1 - y0
  eff <- data.frame(factor = order, effect = unname(effects[order]), stringsAsFactors = FALSE)
  eff$share <- if (delta != 0) eff$effect / delta * 100 else NA_real_
  if (!is.null(formulas)) eff$formula <- unname(formulas[order])
  structure(
    list(
      method = method, model = m, order = order,
      base = d$base, actual = d$actual, n_items = nrow(d$base),
      y0 = y0, y1 = y1, delta = delta,
      effects = eff, residual = delta - sum(eff$effect),
      details = details
    ),
    class = "dfa"
  )
}

.chain_steps <- function(m, d, order) {
  k <- length(order)
  vals <- d$base
  y <- .total(m, vals)
  for (v in order) {
    vals[[v]] <- d$actual[[v]]
    y <- c(y, .total(m, vals))
  }
  labels <- c(
    paste0(m$lhs, "0"),
    if (k > 1L) paste0(m$lhs, "_усл", seq_len(k - 1L)),
    paste0(m$lhs, "1")
  )
  used_actual <- outer(0:k, seq_len(k), ">=")
  tbl <- data.frame(`Расчёт` = labels, check.names = FALSE, stringsAsFactors = FALSE)
  for (j in seq_len(k)) {
    v <- order[j]
    tbl[[v]] <- if (nrow(d$base) == 1L) {
      ifelse(used_actual[, j], d$actual[[v]], d$base[[v]])
    } else {
      ifelse(used_actual[, j], "1", "0")
    }
  }
  tbl[[m$lhs]] <- y
  tbl[["Влияние"]] <- c(NA, diff(y))
  tbl[["Фактор"]] <- c("—", order)
  list(y = y, table = tbl)
}

.gauss_legendre <- function(n) {
  i <- seq_len(n - 1L)
  b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1L)] <- b
  J[cbind(i + 1L, i)] <- b
  e <- eigen(J, symmetric = TRUE)
  list(t = (e$values + 1) / 2, w = e$vectors[1L, ]^2)
}

.fmt <- function(x, digits = 2) {
  out <- formatC(x, format = "f", digits = digits, big.mark = " ")
  out[is.na(x)] <- ""
  out
}

.fmt_signed <- function(x, digits = 2) {
  out <- .fmt(x, digits)
  pos <- !is.na(x) & x > 0
  out[pos] <- paste0("+", out[pos])
  out
}
