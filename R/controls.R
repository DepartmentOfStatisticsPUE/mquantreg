

glmqmControl <- function(maxit = 20, acc = 1e-04, k = 1.345, psi = psi.huber, 
                         sparse = FALSE, eps = .Machine$double.eps^(2/3), lm = 'rlm') 
{
  if (!is.numeric(maxit) || maxit <= 0) 
    stop("maximum number of iterations must be > 0")
  if (!is.numeric(acc) || acc <= 0) 
    stop("value of acc must be > 0")
  if (!is.numeric(k) || k <= 0) 
    stop("value of k must be > 0")
  if (!is.numeric(eps) || eps <= 0) 
    stop("value of eps must be > 0")
  if (sparse) 
    stop("Currently only dense matrices are implemented")
  if (lm != 'rlm') 
    stop("Currently only rlm is supported")
  list(maxit = maxit, acc = acc, k = k, psi = psi, 
       sparse = sparse, eps = eps, lm = lm)
}


rlmControl <- function(init = "ls", scale.est = 'MAD', 
                       method = 'M', test.vec = 'resid') 
{
  if (!init %in% c("ls","lts")) 
    stop("init should equal to ls or lts. See MASS::rlm documentation.")
  if (!scale.est %in% c("MAD", "Huber", "proposal 2")) 
    stop("scale.est should be equal o MAD, Huber or proposal 2. See MASS::rlm documentation.")
  if (!method %in% c("M", "MM")) 
    stop("method should be equal to M or MM. See MASS::rlm documentation.")
  if (test.vec != 'resid') 
    stop("test.vec supports only resid. See MASS::rlm documentation.")
  list(init = init, scale.est = scale.est,
       method = method, test.vec = test.vec)
}


glmrobControl <- function(method = 'Mqle', method.control = glmrobMqle.control(),
                          weights.on.x = 'none', model = FALSE, x = FALSE, y = FALSE,
                          trace.lev = 0, start = NULL) 
{
  if (!method %in% c("Mqle", "BY", "WBY", "MT")) 
    stop("method should be Mqle, BY, WBY or MT. See robustbase::glmrob documentation.")
  if (!is.list(method.control)) 
    stop("method.control should be either glmrobMqle.control, glmrobBY.control 
         or glmrobMT.control. See robustbase::glmrob..control documentation.")
  if (!weights.on.x %in% c("none", "hat", "robCov", "covMcd")) 
    stop("method should be none, hat, robCov or covMcd. See robustbase::glmrob documentation.")
  if (!is.logical(model) || !is.logical(x) || !is.logical(y)) 
    stop("Arguments model, x and y should be logical. See robustbase::glmrob documentation.")
  if (!is.numeric(trace.lev) || trace.lev < 0) 
    stop("trace.lev must be >= 0. See robustbase::glmrob documentation.")
  list(method = method, method.control = method.control,
       weights.on.x = weights.on.x, model = model, x = x, y = y,
       trace.lev = trace.lev, start = start)
}


