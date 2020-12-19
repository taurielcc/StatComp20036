#include <Rcpp.h>
#include <cmath>
#include <unordered_map>
#include <cstring>
using namespace Rcpp;

NumericMatrix L2SqrDistance(NumericMatrix X, NumericMatrix Z)
{
    auto n = X.nrow();
    auto m = Z.nrow();
    auto d = X.ncol();

    NumericMatrix result(m, n);

    for (size_t i = 0; i < m; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            double s = 0;
            for (size_t k = 0; k < d; ++k)
            {
                s += (X(j, k) - Z(i, k)) * (X(j, k) - Z(i, k));
            }
            result(i, j) = s;
        }
    }

    return result;
}

template <typename T>
struct indexed_data_t
{
    T data;
    size_t index;

    bool operator<(const indexed_data_t &other)
    {
        return this->data < other.data;
    }
};

template <typename T>
void heapify(T *arr, int idx, size_t n)
{
    while (idx < n)
    {
        auto left_idx = idx * 2 + 1;
        auto right_idx = idx * 2 + 2;
        auto target_idx = idx;

        if (left_idx >= n)
            break;
        if (arr[left_idx] < arr[idx] && (right_idx >= n || arr[left_idx] < arr[right_idx]))
        {
            target_idx = left_idx;
        }
        else if (right_idx < n && arr[right_idx] < arr[idx])
        {
            target_idx = right_idx;
        }
        if (target_idx != idx)
        {
            std::swap(arr[idx], arr[target_idx]);
            idx = target_idx;
        }
        else
        {
            break;
        }
    }
}

IntegerMatrix topK(NumericMatrix dist, int k)
{
    auto n_samples = dist.nrow();
    auto n_dim = dist.ncol();
    auto temp_arr = new indexed_data_t<double>[n_dim];

    IntegerMatrix result(n_samples, k);

    for (size_t i = 0; i < n_samples; i++)
    {
        for (size_t j = 0; j < n_dim; j++)
        {
            temp_arr[j].index = j;
            temp_arr[j].data = dist(i, j);
        }

        for (int j = n_dim / 2 - 1; j >= 0; j--)
        {
            heapify(temp_arr, j, n_dim);
        }

        for (int j = 0; j < k; j++)
        {
            result(i, j) = temp_arr[0].index;
            std::swap(temp_arr[0], temp_arr[n_dim - j - 1]);
            heapify(temp_arr, 0, n_dim - j - 1);
        }
    }

    delete[] temp_arr;
    return result;
}

IntegerVector selectTopK(IntegerMatrix idx_mat, IntegerVector labels, int n_classes)
{
    auto n = idx_mat.nrow();
    auto k = idx_mat.ncol();
    auto cnt = new int[n_classes];

    IntegerVector result(n);
    for (int i = 0; i < n; ++i)
    {
        std::memset(cnt, 0, sizeof(int) * n_classes);
        for (int j = 0; j < k; ++j)
            cnt[labels[idx_mat(i, j)]]++;

        int most_freq_label = -1, max_freq = 0;
        for (int j = 0; j < n; ++j)
            if (max_freq < cnt[j])
            {
                max_freq = cnt[j];
                most_freq_label = j;
            }

        result(i) = most_freq_label;
    }

    delete[] cnt;
    return result;
}


//' @title A KNN function using Rcpp
//' @description A KNN(K-Nearest Neighbor) classification for test set from training set using Rcpp. For each row of the test set, the k nearest (in Euclidean distance) training set vectors are found, and the classification is decided by majority vote.
//' @param X matrix of training set cases
//' @param y label of true classification of training set
//' @param Z matrix of test set cases
//' @param n_classes number of classes in y
//' @param k number of neighbours considerd
//' @return dist Euclidean distance matrix
//' @return results label of classification of test set
//' @examples
//' \dontrun{
//' X<-matrix(c(1:12),ncol=2)
//' y<-c(0,0,0,1,1,1)
//' Z<-matrix(c(2,6,7,11),ncol=2)
//' res<-KNN(X,y,Z,2,k=3)
//' }
//' @export
// [[Rcpp::export]]
List kNN(NumericMatrix X, IntegerVector y, NumericMatrix Z, int n_classes, int k)
{
    auto dist = L2SqrDistance(X, Z);
    auto topk_idx = topK(dist, k);
    auto pred_y = selectTopK(topk_idx, y, n_classes);
    return List::create(Named("dist") = dist, Named("results") = pred_y);
}