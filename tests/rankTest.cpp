#include<iostream>
#include"HDF5Reader/HDF5Reader.h"
#include"rank/rankNormalize.h"
using namespace std;

int main() {
	HDF5reader hr;
	hr.readHDF5File("C:/Users/ltw/Desktop/ica_cord_blood_h5.h5");
	rankNormalize rn(hr);
	rn.ranks();
	rn.print();
}
