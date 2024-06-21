// [[Rcpp::depends(TMB)]]
#include "RTMB.h"

// Workarounds needed for parallel case:
// R_CheckStack() is junk when called from non-master thread.
// We can effectively disable the check by setting 'R_CStackDir = 0'.
#ifdef _OPENMP
extern int R_CStackDir;
#endif
struct CStackWorkaround {
#ifdef _OPENMP
  int old;
  void begin() {
    old = R_CStackDir;
    if (omp_get_thread_num() != 0) R_CStackDir = 0;
  }
  void end() {
    R_CStackDir = old;
  }
#else
  void begin() { }
  void end() { }
#endif
};

namespace TMBad {
/*
  EvalOp: Operator that evaluates an R function taking single numeric
  scalar as input.

  Memory management:

  - 'Rcpp::Function' automatically PROTECTs/UNPROTECTs. This is in
    general useful for temporary liftime objects, but we have to be
    careful when placing a 'Rcpp::Function' as a member of an object
    with 'unlimited lifetime' (e.g. an AD operator). If the object is
    copied (from R), it'll invoke PROTECTs that are not followed by
    the corresponding UNPROTECTs. A simple solution is to place
    'Rcpp::Function' in a shared pointer. That effectively disables
    the Rcpp memory management, while still taking advantage of the
    Rcpp cleanup code triggered when the last reference dies. I have
    verified that cleanup code is indeed triggered when the tapes
    (obj) *and* R function (F) are removed from the R
    workspace. Removing just F is not enough (nice).

  - Evaluating F from C++ is known to not work in parallel. We use
    'omp critical' to make sure F is never evaluated by two threads at
    the same time. However, that's not enough. R occasionally runs
    'R_StackCheck()' which works in a highly thead unsafe manner. It
    will fail for any other thread than the master. Inspection of
    'R_StackCheck()' reveals that it can be disabled by temporarily
    setting 'R_CStackDir = 0'. FIXME: We could probably work out the
    rest of these 'R_CStack...' variables in order for
    'R_StackCheck()' to work reliably for all threads. Don't know if
    this is worth while though.
*/
template <bool with_derivs = false>
struct EvalOp : global::DynamicOperator< -1 , -1 > {
  static const bool have_input_size_output_size = true;
  static const bool add_forward_replay_copy = true;
  std::shared_ptr<Rcpp::Function> Fptr; // forward
  std::shared_ptr<Rcpp::Function> Rptr; // reverse
  size_t m, n;
  Index input_size()  const { return m; }
  Index output_size() const { return n; }
  EvalOp (Rcpp::Function F, size_t m, size_t n) :
    Fptr(std::make_shared<Rcpp::Function>(F)),
    m(m),
    n(n) {
    if (with_derivs) {
      Rptr = std::make_shared<Rcpp::Function>(F.attr("reverse"));
    }
  }
  void forward(ForwardArgs<double> &args) {
#ifdef _OPENMP
#pragma omp critical
    {
#endif
      CStackWorkaround R;
      R.begin();
      Rcpp::NumericVector i(m);
      for (size_t l=0; l<m; l++) i[l] = args.x(l);
      SEXP y = (*Fptr)(i);
      // FIXME: Any Rcpp way of handling arbitrary output? For now doing PROTECT/UNPROTECT manually...
      PROTECT(y);
      if ((size_t) LENGTH(y) != n) {
        R.end();
	UNPROTECT(1); // y
        Rcpp::stop("Wrong output length");
      }
      if (Rf_isReal(y)) {
        double* py = REAL(y);
        for (size_t i=0; i<n; i++) { args.y(i) = py[i]; }
      } else if (Rf_isInteger(y)) {
        int* py = INTEGER(y);
        for (size_t i=0; i<n; i++) { args.y(i) = py[i]; }
      } else {
        R.end();
	UNPROTECT(1); // y
        Rcpp::stop("EvalOp: Function must return 'real' or 'integer'");
      }
      R.end();
      UNPROTECT(1); // y
#ifdef _OPENMP
    }
#endif
  }
  template <class Type> void forward(ForwardArgs<Type> &args) {
    TMBAD_ASSERT(false);
  }
  void reverse(ReverseArgs<double> &args) {
    if (with_derivs) {
      Rcpp::NumericVector x(m);
      Rcpp::NumericVector y(n);
      Rcpp::NumericVector dy(n);
      for (size_t l=0; l<m; l++) {
        x[l] = args.x(l);
      }
      for (size_t l=0; l<n; l++) {
        y[l] = args.y(l);
        dy[l] = args.dy(l);
      }
      Rcpp::NumericVector wtJ = (*Rptr)(x,y,dy);
      if ( (size_t) wtJ.size() != m)
        Rcpp::stop("Wrong length of 'reverse(x,y,dy)' = t(dy) %*% jacobian(x)");
      for (size_t l=0; l<m; l++) args.dx(l) += wtJ[l];
    }
    // otherwise void derivs
  }
  void reverse(ReverseArgs<ad> &args) {
    if (with_derivs) {
      // Exceptionally we add Rcpp try/catch blocks here to prevent
      // crash when called from TMB::MakeADFun (which does not use
      // Rcpp)
      BEGIN_RCPP
      Rcpp::ComplexVector x(m); x = as_advector(x);
      Rcpp::ComplexVector y(n); y = as_advector(y);
      Rcpp::ComplexVector dy(n); dy = as_advector(dy);
      for (size_t l=0; l<m; l++) {
        x[l] = ad2cplx(args.x(l));
      }
      for (size_t l=0; l<n; l++) {
        y[l] = ad2cplx(args.y(l));
        dy[l] = ad2cplx(args.dy(l));
      }
      Rcpp::ComplexVector wtJ = (*Rptr)(x, y, dy); // User code could throw !
      CHECK_INPUT(wtJ); // Check result from R
      if ( (size_t) wtJ.size() != m)
        Rcpp::stop("'%s': Length of derivative (%u) not as expected (%u)",
                   op_name(),
                   (size_t) wtJ.size(),
                   (size_t) m);
      for (size_t l=0; l<m; l++) args.dx(l) += cplx2ad(wtJ[l]);
      VOID_END_RCPP
    }
  }
  template <class Type> void reverse(ReverseArgs<Type> &args) {
    TMBAD_ASSERT(false);
  }
  const char* op_name() {
    SEXP name = Fptr -> attr("name");
    if (name != R_NilValue)
      return CHAR(STRING_ELT(name, 0));
    else
      return "EvalOp";
  }
  void print(TMBad::global::print_config cfg) {
    Rcout << cfg.prefix;
    Rcout << "F=" << *Fptr << " ";
    Rcout << "n=" << n << "\n";
  }
};
}

// [[Rcpp::export]]
Rcpp::ComplexVector TapedEval(Rcpp::Function F, Rcpp::ComplexVector i) {
  if (!ad_context()) Rcpp::stop("TapedSubset requires an active ad context");
  CHECK_INPUT(i);
  size_t m = i.size();
  ad* pi = adptr(i);
  // Test eval to get n
  Rcpp::NumericVector i_test(m);
  for (size_t l=0; l<m; l++) i_test[l] = pi[l].Value();
  Rcpp::NumericVector y_test = F(i_test);
  size_t n = LENGTH(y_test);
  // Add to tape
  std::vector<ad> x(pi, pi + m);
  bool with_derivs = F.hasAttribute("reverse");
  std::vector<ad> y = (with_derivs ?
                       TMBad::global::Complete<TMBad::EvalOp<true > >(F, m, n) (x) :
                       TMBad::global::Complete<TMBad::EvalOp<false> >(F, m, n) (x) );
  // Pass to R
  Rcpp::ComplexVector ans(n);
  for (size_t j=0; j < n; j++) {
    ans[j] = ad2cplx(y[j]);
  }
  DUPLICATE_ATTRIB(ans, y_test);
  return as_advector(ans);
}

/* Interface to some computational graph transforms which are also available from TMB*/

// Helper for 'laplace' and 'newton'
void remove_random_parameters(TMBad::ADFun<>* adf, const std::vector<TMBad::Index>& random) {
  std::vector<bool> mask(adf->Domain(), true);
  for (size_t i = 0; i<random.size(); i++)
    mask[random[i]] = false;
  adf->glob.inv_index = TMBad::subset(adf->glob.inv_index, mask);
}
std::vector<TMBad::Index> zero_based_unique_index (const std::vector<TMBad::Index> &x, TMBad::Index max) {
  std::vector<TMBad::Index> y(x);
  std::vector<bool> mark(max, false);
  for (size_t i=0; i<y.size(); i++) {
    y[i]--;
    if (y[i] >= max) Rcpp::stop("Index out of bounds");
    if (mark[y[i]])  Rcpp::stop("Index not unique");
    mark[y[i]] = true;
  }
  return y;
}
void laplace_transform(TMBad::ADFun<>* adf, std::vector<TMBad::Index> random, SEXP config) {
  if (random.size() == 0) return;
  random = zero_based_unique_index(random, adf->Domain());
  newton::newton_config cfg(config);
  *adf = newton::Laplace_(*adf, random, cfg);
  remove_random_parameters(adf, random);
}
void newton_transform(TMBad::ADFun<>* adf, std::vector<TMBad::Index> random, SEXP config) {
  // TMB FIXME: This code is almost copy-paste from 'newton::Laplace_'. Add it 'Newton_' there.
  if (random.size() == 0) return;
  random = zero_based_unique_index(random, adf->Domain());;
  newton::newton_config cfg(config);
  newton::slice<> S(*adf, random);
  TMBad::ADFun<> ans;
  std::vector<double> xd = (*adf).DomainVec();
  S.x = std::vector<ad> (xd.begin(), xd.end());
  ans.glob.ad_start();
  TMBad::Independent(S.x);
  vector<ad> start = TMBad::subset(S.x, random);
  std::vector<ad> y = newton::Newton(S, start, cfg);
  TMBad::Dependent(y);
  ans.glob.ad_stop();
  *adf = ans;
  remove_random_parameters(adf, random);
}

// Set low-rank tags
// [[Rcpp::export]]
Rcpp::ComplexVector LowRankTag(const Rcpp::ComplexVector &x) {
  CHECK_INPUT(x);
  size_t n = x.size();
  Rcpp::ComplexVector y(n);
  ad* X = adptr(x);
  ad* Y = adptr(y);
  for (size_t i=0; i<n; i++) Y[i] = newton::Tag(X[i]);
  return as_advector(y);
}
