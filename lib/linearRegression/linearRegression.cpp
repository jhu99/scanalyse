#include"linearRegression.h"

template<class T> LinearRegression<T>::LinearRegression() {

}

template<class T> LinearRegression<T>::LinearRegression(Cells<T> cells) {
	this->n = cells.getGeneToNum().size();
	this->p = cells.getCellToNum().size();
	this->cells = cells;
	linearRegressionParameters = new LinearRegressionParameter *[p];
	for (int i = 0; i < p; i++) {
		linearRegressionParameters[i] = new LinearRegressionParameter[n];
	}
}


template <class T> LinearRegression<T>::~LinearRegression() {
	for (int i = 0; i < p; i++) {
		delete[] linearRegressionParameters[i];
	}
	delete[] linearRegressionParameters;
}

template <class T> void LinearRegression<T>::setCell(Cells<T> cells) {
	this->cells = cells;
}
template <class T> void LinearRegression<T>::setLinearRegressionParameter(LinearRegressionParameter **linearRegressionParameters) {
	this->linearRegressionParameters = linearRegressionParameters;
}

template <class T> Cells<T> LinearRegression<T>::getCell() {
	return cells;
}
template <class T> LinearRegressionParameter **LinearRegression<T>::getLinearRegressionParameter() {
	return linearRegressionParameters;
}
template <class T>void LinearRegression<T>::calculate(int i) {
	double c0 = 0, c1 = 0, cov00 = 0, cov01 = 0, cov11 = 0, sumsq = 0, sumtot = 0;

	for (int j = 0; j < p; j++) {
		if(j==i)continue;
		mut.lock();
		clock_t start = clock();
		double *x, *y;
		x = new double[n];
		y = new double[n];
		//memcpy(x, cells.getCell()[i], sizeof(x));
		//memcpy(y, cells.getCell()[j], sizeof(y));
		for (int k = 0; k < n; k++) {
			x[k] = cells.getCell()[i][k];
			y[k] = cells.getCell()[j][k];
		}

		//sumtot = gsl_stats_tss(y, 1, n);
		sumtot =  gsl_stats_variance(y,1,n)*(n-1);
		gsl_fit_linear(x, 1, y, 1, n, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
		
		LinearRegressionParameter lrp(c0, c1, cov00, cov01, cov11, sumsq,sumtot,n);
		//lrp.print();
		linearRegressionParameters[i][j] = lrp;

		delete[] x;
		delete[] y;
		clock_t ends = clock();
		cout <<"Running Time : "<<(double)(ends - start)/ CLOCKS_PER_SEC << endl;
		mut.unlock();
	}
}
template <class T> void LinearRegression<T>::calculateLinearRegression(int nt) {
	thread *threads;
	threads = new thread[nt];
	for (int k = 0; k < p / nt + 1; k++) {
		for (int i = 0; i < min(nt, p - nt * k); i++) {
			//cout << k * nt + i << endl;
			threads[i] = thread(&LinearRegression<T>::calculate, this, k*nt + i);
		}

		for (int i = 0; i < min(nt, p - nt * k); i++) {
			threads[i].join();
		}
	}
	delete[] threads;
}

template <class T> void LinearRegression<T>::writeFile(string path) {
	ofstream oFile;

	oFile.open(path, ios::out | ios::trunc);

	unordered_map<int, string> cls = cells.getNumToCell();
	oFile << ",";
	for (int i = 0; i < p; i++) {
		oFile << cls[i] << ",";
	}
	oFile << endl;
	for (int i = 0; i < p; i++) {
		oFile << cls[i] << ",";
		for (int j = 0; j < p; j++) {
			if (i == j) {
				oFile << ",";
				continue;
			}
			string s = to_string(linearRegressionParameters[i][j].getC0()) + "|" + to_string(linearRegressionParameters[i][j].getC1()) 
				+ "|" + to_string(linearRegressionParameters[i][j].getCov00()) + "|"+ to_string(linearRegressionParameters[i][j].getCov01()) 
				+ "|" + to_string(linearRegressionParameters[i][j].getCov11()) + "|" + to_string(linearRegressionParameters[i][j].getSumsq()) 
				+ "|"+ to_string(linearRegressionParameters[i][j].getSumtot()) + "|" + to_string(linearRegressionParameters[i][j].getR2()) 
				+ "|" + to_string(linearRegressionParameters[i][j].getT()) + "|"+ to_string(linearRegressionParameters[i][j].getPvalue());
			oFile << s << ",";
		}
		oFile << endl;
	}
	oFile.close();
}
