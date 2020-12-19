## -----------------------------------------------------------------------------
boot <- function(data,func=NULL, B){
  theta.hat <- func(data)
  #set up the bootstrap
  n <- length(data)      #sample size
  theta.b <- numeric(B)     #storage for replicates
  for (b in 1:B) {
    #randomly select the indices
    i <- sample(1:n, size = n, replace = TRUE)
    dat <- data[i]       #i is a vector of indices
    theta.b[b] <- func(dat)
  }
  #bootstrap estimate of standard error of R
  bias.theta <- mean(theta.b - theta.hat)
  se <- sd(theta.b)
  return(list(bias.b = bias.theta,se.b = se))
}

## -----------------------------------------------------------------------------
jack <- function(data,func=NULL){
  theta.hat <- func(data)
  #set up the bootstrap
  #B is the number of replicates
  n <- length(data)      #sample size
  M <- numeric(n)
  for (i in 1:n) { #leave one out
    y <- data[-i]
    M[i] <- func(y)
  }
  Mbar <- mean(M)
  se.jack <- sqrt(((n - 1)/n) * sum((M - Mbar)^2))
  return(se.jack)
}

## -----------------------------------------------------------------------------
jackafterboot <- function(data,func=NULL,B){
  n <- length(data)
  theta.b <- numeric(B)
  # set up storage for the sampled indices
  indices <- matrix(0, nrow = B, ncol = n)
  # jackknife-after-bootstrap step 1: run the bootstrap
  for (b in 1:B) {
    i <- sample(1:n, size = n, replace = TRUE)
    y <- data[i]
    theta.b[b] <- func(y)
    #save the indices for the jackknife
    indices[b, ] <- i
  }
  #jackknife-after-bootstrap to est. se(se)
  se.jack <- numeric(n)
  for (i in 1:n) {
    #in i-th replicate omit all samples with x[i]
    keep <- (1:B)[apply(indices, MARGIN = 1,
                        FUN = function(k) {!any(k == i)})]
    se.jack[i] <- sd(theta.b[keep])
  }
  se.boot <- sd(theta.b)
  se.jackafterboot <- sqrt((n-1) * mean((se.jack - mean(se.jack))^2))
  return(list(se.boot = se.boot, se.jackafterboot=se.jackafterboot))
}

## -----------------------------------------------------------------------------
data <- 20 * rbeta(1000,2,3)
(boot <- boot(data = data, B = 200, func = median)$se.b)
(jack <- jack(data = data, func = median))
(jackafterboot<- jackafterboot(data, B = 200, func = median))

## ----eval=FALSE---------------------------------------------------------------
#  NumericMatrix L2SqrDistance(NumericMatrix X, NumericMatrix Z)
#  {
#      auto n = X.nrow();
#      auto m = Z.nrow();
#      auto d = X.ncol();
#  
#      NumericMatrix result(m, n);
#  
#      for (size_t i = 0; i < m; ++i)
#      {
#          for (size_t j = 0; j < n; ++j)
#          {
#              double s = 0;
#              for (size_t k = 0; k < d; ++k)
#              {
#                  s += (X(j, k) - Z(i, k)) * (X(j, k) - Z(i, k));
#              }
#              result(i, j) = s;
#          }
#      }
#  
#      return result;
#  }
#  
#  struct indexed_data_t
#  {
#      T data;
#      size_t index;
#  
#      bool operator<(const indexed_data_t &other)
#      {
#          return this->data < other.data;
#      }
#  };
#  
#  void heapify(T *arr, int idx, size_t n)
#  {
#      while (idx < n)
#      {
#          auto left_idx = idx * 2 + 1;
#          auto right_idx = idx * 2 + 2;
#          auto target_idx = idx;
#  
#          if (left_idx >= n)
#              break;
#          if (arr[left_idx] < arr[idx] && (right_idx >= n || arr[left_idx] < arr[right_idx]))
#          {
#              target_idx = left_idx;
#          }
#          else if (right_idx < n && arr[right_idx] < arr[idx])
#          {
#              target_idx = right_idx;
#          }
#          if (target_idx != idx)
#          {
#              std::swap(arr[idx], arr[target_idx]);
#              idx = target_idx;
#          }
#          else
#          {
#              break;
#          }
#      }
#  }
#  
#  IntegerMatrix topK(NumericMatrix dist, int k)
#  {
#      auto n_samples = dist.nrow();
#      auto n_dim = dist.ncol();
#      auto temp_arr = new indexed_data_t<double>[n_dim];
#  
#      IntegerMatrix result(n_samples, k);
#  
#      for (size_t i = 0; i < n_samples; i++)
#      {
#          for (size_t j = 0; j < n_dim; j++)
#          {
#              temp_arr[j].index = j;
#              temp_arr[j].data = dist(i, j);
#          }
#  
#          for (int j = n_dim / 2 - 1; j >= 0; j--)
#          {
#              heapify(temp_arr, j, n_dim);
#          }
#  
#          for (int j = 0; j < k; j++)
#          {
#              result(i, j) = temp_arr[0].index;
#              std::swap(temp_arr[0], temp_arr[n_dim - j - 1]);
#              heapify(temp_arr, 0, n_dim - j - 1);
#          }
#      }
#  
#      delete[] temp_arr;
#      return result;
#  }
#  
#  IntegerVector selectTopK(IntegerMatrix idx_mat, IntegerVector labels, int n_classes)
#  {
#      auto n = idx_mat.nrow();
#      auto k = idx_mat.ncol();
#      auto cnt = new int[n_classes];
#  
#      IntegerVector result(n);
#      for (int i = 0; i < n; ++i)
#      {
#          std::memset(cnt, 0, sizeof(int) * n_classes);
#          for (int j = 0; j < k; ++j)
#              cnt[labels[idx_mat(i, j)]]++;
#  
#          int most_freq_label = -1, max_freq = 0;
#          for (int j = 0; j < n; ++j)
#              if (max_freq < cnt[j])
#              {
#                  max_freq = cnt[j];
#                  most_freq_label = j;
#              }
#  
#          result(i) = most_freq_label;
#      }
#  
#      delete[] cnt;
#      return result;
#  }
#  
#  List kNN(NumericMatrix X, IntegerVector y, NumericMatrix Z, int n_classes, int k)
#  {
#      auto dist = L2SqrDistance(X, Z);
#      auto topk_idx = topK(dist, k);
#      auto pred_y = selectTopK(topk_idx, y, n_classes);
#      return List::create(Named("dist") = dist, Named("results") = pred_y);
#  }

