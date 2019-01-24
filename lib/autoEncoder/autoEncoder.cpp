#include"autoEncoder.h"

AutoEncoder::AutoEncoder(int n_input,int size){
	this->n_input = n_input;
	this->n_output = n_input;
	this->size = size;
	this->n_hidden1 = 128;
	this->n_hidden2 = 32;
	this->n_hidden3 = 2;
	this->n_hidden4 = n_hidden2;
	this->n_hidden5 = n_hidden1;
}

AutoEncoder::~AutoEncoder(){

}

vector<float> AutoEncoder::getData(double** data,int size){
	vector<float> v;
	for(int i=0;i<size;i++){
		for(int j=0;j<n_input;j++){
			v.push_back(data[i][j]);
		}	
	}
	return v;
}

void AutoEncoder::train(int n_train,string path_read, string type){
			
	Tensor x_datas(DataTypeToEnum<float>::v(),TensorShape{size, n_input});
  	//copy_n(x_data.begin(), x_data.size(),x_datas.flat<float>().data());
	Tensor y_datas(DataTypeToEnum<float>::v(),TensorShape{size, n_input});
  	//copy_n(y_data.begin(), y_data.size(),y_datas.flat<float>().data());
  
  	Scope scope = Scope::NewRootScope();

  	auto x = Placeholder(scope, DT_FLOAT);
  	auto y = Placeholder(scope, DT_FLOAT);

  	// weights init
  	auto w1 = Variable(scope, {n_input, n_hidden1}, DT_FLOAT);
  	auto assign_w1 = Assign(scope, w1, RandomNormal(scope, {n_input, n_hidden1}, DT_FLOAT));

  	auto w2 = Variable(scope, {n_hidden1, n_hidden2}, DT_FLOAT);
  	auto assign_w2 = Assign(scope, w2, RandomNormal(scope, {n_hidden1, n_hidden2}, DT_FLOAT));

	auto w3 = Variable(scope, {n_hidden2, n_hidden3}, DT_FLOAT);
  	auto assign_w3 = Assign(scope, w3, RandomNormal(scope, {n_hidden2, n_hidden3}, DT_FLOAT));

	auto w4 = Variable(scope, {n_hidden3, n_hidden4}, DT_FLOAT);
  	auto assign_w4 = Assign(scope, w4, RandomNormal(scope, {n_hidden3, n_hidden4}, DT_FLOAT));

	auto w5 = Variable(scope, {n_hidden4, n_hidden5}, DT_FLOAT);
  	auto assign_w5 = Assign(scope, w5, RandomNormal(scope, {n_hidden4, n_hidden5}, DT_FLOAT));

	auto w6 = Variable(scope, {n_hidden5, n_output}, DT_FLOAT);
  	auto assign_w6 = Assign(scope, w6, RandomNormal(scope, {n_hidden5, n_output}, DT_FLOAT));


  	// bias init
  	auto b1 = Variable(scope, {1, n_hidden1}, DT_FLOAT);
  	auto assign_b1 = Assign(scope, b1, RandomNormal(scope, {1,  n_hidden1}, DT_FLOAT));

  	auto b2 = Variable(scope, {1,n_hidden2}, DT_FLOAT);
  	auto assign_b2 = Assign(scope, b2, RandomNormal(scope, {1, n_hidden2}, DT_FLOAT));
			
	auto b3 = Variable(scope, {1, n_hidden3}, DT_FLOAT);
  	auto assign_b3 = Assign(scope, b3, RandomNormal(scope, {1, n_hidden3}, DT_FLOAT));

	auto b4 = Variable(scope, {1, n_hidden4}, DT_FLOAT);
  	auto assign_b4 = Assign(scope, b4, RandomNormal(scope, {1, n_hidden4}, DT_FLOAT));

	auto b5 = Variable(scope, {1, n_hidden5}, DT_FLOAT);
  	auto assign_b5 = Assign(scope, b5, RandomNormal(scope, {1, n_hidden5}, DT_FLOAT));

	auto b6 = Variable(scope, {1, n_output}, DT_FLOAT);
  	auto assign_b6 = Assign(scope, b6, RandomNormal(scope, {1, n_output}, DT_FLOAT));

  	// layers
  	auto layer_1 = Tanh(scope, Add(scope, MatMul(scope, x, w1), b1));
  	auto layer_2 = Tanh(scope, Add(scope, MatMul(scope, layer_1, w2), b2));
	auto layer_3 = Tanh(scope, Add(scope, MatMul(scope, layer_2, w3), b3));
  	auto layer_4 = Tanh(scope, Add(scope, MatMul(scope, layer_3, w4), b4));
	auto layer_5 = Tanh(scope, Add(scope, MatMul(scope, layer_4, w5), b5));
  	auto layer_6 = Tanh(scope, Add(scope, MatMul(scope, layer_5, w6), b6));

  	// regularization
  	auto regularization = AddN(scope,
                        initializer_list<Input>{L2Loss(scope, w1),
                                                L2Loss(scope, w2),
						L2Loss(scope, w3),
						L2Loss(scope, w4),
						L2Loss(scope, w5),
						L2Loss(scope, w6)});

  	// loss calculation
  	auto loss = Add(scope,
                  ReduceMean(scope, Square(scope, Sub(scope, layer_6, y)), {0,1}),
                  Mul(scope, Cast(scope, 0.01,  DT_FLOAT), regularization));

  	// add the gradients operations to the graph
  	std::vector<Output> grad_outputs;
  	TF_CHECK_OK(AddSymbolicGradients(scope, {loss}, {w1,w2,w3,w4,w5,w6, b1,b2,b3,b4,b5,b6}, &grad_outputs));

  	// update the weights and bias using gradient descent
  	auto apply_w1 = ApplyGradientDescent(scope, w1, Cast(scope, 0.01,  DT_FLOAT), {grad_outputs[0]});
  	auto apply_w2 = ApplyGradientDescent(scope, w2, Cast(scope, 0.01,  DT_FLOAT), {grad_outputs[1]});
	auto apply_w3 = ApplyGradientDescent(scope, w3, Cast(scope, 0.01,  DT_FLOAT), {grad_outputs[2]});
  	auto apply_w4 = ApplyGradientDescent(scope, w4, Cast(scope, 0.01,  DT_FLOAT), {grad_outputs[3]});
	auto apply_w5 = ApplyGradientDescent(scope, w5, Cast(scope, 0.01,  DT_FLOAT), {grad_outputs[4]});
  	auto apply_w6 = ApplyGradientDescent(scope, w6, Cast(scope, 0.01,  DT_FLOAT), {grad_outputs[5]});

  	auto apply_b1 = ApplyGradientDescent(scope, b1, Cast(scope, 0.01,  DT_FLOAT), {grad_outputs[6]});
  	auto apply_b2 = ApplyGradientDescent(scope, b2, Cast(scope, 0.01,  DT_FLOAT), {grad_outputs[7]});
	auto apply_b3 = ApplyGradientDescent(scope, b3, Cast(scope, 0.01,  DT_FLOAT), {grad_outputs[8]});
  	auto apply_b4 = ApplyGradientDescent(scope, b4, Cast(scope, 0.01,  DT_FLOAT), {grad_outputs[9]});
	auto apply_b5 = ApplyGradientDescent(scope, b5, Cast(scope, 0.01,  DT_FLOAT), {grad_outputs[10]});
  	auto apply_b6 = ApplyGradientDescent(scope, b6, Cast(scope, 0.01,  DT_FLOAT), {grad_outputs[11]});

			

  	ClientSession session(scope);
  	std::vector<Tensor> outputs;
  
  	// init the weights and biases by running the assigns nodes once
  	TF_CHECK_OK(session.Run({assign_w1, assign_w2,assign_w3, assign_w4,assign_w5, assign_w6,assign_b1, assign_b2,assign_b3, assign_b4,assign_b5, assign_b6}, nullptr));
  
  	// training steps
	SparseMatrix sm;
	sm.readHDF5File(path_read, type);
  	for (int i = 0; i < n_train; ++i) {
		x_data = getData(sm.fetch_batch(i,size),size);
		y_data = x_data;
		copy_n(x_data.begin(), x_data.size(),x_datas.flat<float>().data());
  		copy_n(y_data.begin(), y_data.size(),y_datas.flat<float>().data());
    		if (i % 2 == 0) {
      			TF_CHECK_OK(session.Run({{x, x_datas}, {y, y_datas}}, {loss},  &outputs));
      			std::cout << "Loss after " << i << " steps " << outputs[0].scalar<float>() << std::endl;
    		}
   	 	// nullptr because the output from the run is useless
    		TF_CHECK_OK(session.Run({{x, x_datas}, {y, y_datas}}, {apply_w1, apply_w2,apply_w3, apply_w4,apply_w5, apply_w6, apply_b1, apply_b2,apply_b3, apply_b4,apply_b5, apply_b6}, nullptr));
	}
} 
