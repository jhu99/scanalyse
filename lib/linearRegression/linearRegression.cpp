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
	double c0 = 0, c1 = 0, cov00 = 0, cov01 = 0, cov11 = 0, sumsq = 0;

	for (int j = 0; j < p; j++) {
		mut.lock();
		double *x, *y;
		x = new double[p];
		y = new double[p];
		memcpy(x, cells.getCell()[i], sizeof(x));
		memcpy(y, cells.getCell()[j], sizeof(y));

		gsl_fit_linear(x, 1, y, 1, p, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);

		LinearRegressionParameter lrp(c0, c1, cov00, cov01, cov11, sumsq);
		linearRegressionParameters[i][j] = lrp;
		//cout << c0 << ": " << c1 << ":" << nn << " " << i << endl;
		nn++;
		delete[] x;
		delete[] y;
		mut.unlock();
	}
}
template <class T> void LinearRegression<T>::calculateLinearRegression(int nt) {
	thread *threads;
	threads = new thread[nt];
	for (int k = 0; k < p / nt + 1; k++) {
		for (int i = 0; i < min(nt, p - nt * k); i++) {
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
			string s = to_string(linearRegressionParameters[i][j].getC0()) + "|" + to_string(linearRegressionParameters[i][j].getC1()) + "|" + to_string(linearRegressionParameters[i][j].getCov00()) + "|" + to_string(linearRegressionParameters[i][j].getCov01()) + "|" + to_string(linearRegressionParameters[i][j].getCov11()) + "|" + to_string(linearRegressionParameters[i][j].getSumsq());
			oFile << s << ",";
		}
		oFile << endl;
	}
	oFile.close();
}
