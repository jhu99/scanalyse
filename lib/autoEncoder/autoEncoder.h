#include "tensorflow/cc/client/client_session.h"
#include "tensorflow/cc/ops/standard_ops.h"
#include "tensorflow/core/framework/tensor.h"
#include "tensorflow/cc/framework/gradients.h"
#include "SparseMatrix/SparseMatrix.h"
#include"rank/rankNormalize.h"
#include"qqNorm/caculateInterface.h"
#include <string>  

using namespace tensorflow;
using namespace tensorflow::ops;
using namespace std;

class AutoEncoder{
	private:
		int n_input,n_output;
		int n_hidden1;
		int n_hidden2;
		int n_hidden3;
		int n_hidden4;
		int n_hidden5;
		int size;
		vector<float> x_data,y_data;
	public:
		AutoEncoder(int n_input,int size);
		~AutoEncoder();
		vector<float> getData(double** data,int size);
		void train(int n_train,string path_read, string type);
};
