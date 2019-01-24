#include"autoEncoder/autoEncoder.h"

using namespace tensorflow;
using namespace tensorflow::ops;
using namespace std;

int main(){
	 
	AutoEncoder ae(38400,128);
	ae.train(1000,"./data/ica_cord_blood_h5.h5","original");
}
