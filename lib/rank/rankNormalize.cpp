#include "rankNormalize.h"

rankNormalize::rankNormalize() {

}
rankNormalize::rankNormalize(HDF5reader hr) {
	this->hr = hr;
	n = hr.get_data_count();
	geneInfos = new geneInfo[n];
	rank = new unsigned short[n];
	setGeneInfoByHr();
}
rankNormalize::~rankNormalize() {
	delete[] geneInfos;
	delete[] rank;
	//hr.deleteHDF5();
}
void rankNormalize::setN(long long n) {
	this->n = n;
}
void rankNormalize::setHr(HDF5reader hr) {
	this->hr = hr;
}
void rankNormalize::setGeneInfos(geneInfo *geneInfos) {
	this->geneInfos = geneInfos;
}
void rankNormalize::setRank(unsigned short *rank) {
	this->rank = rank;
}

long long rankNormalize::getN() {
	return n;
}
HDF5reader rankNormalize::getHr() {
	return hr;
}
geneInfo* rankNormalize::getGeneInfos() {
	return geneInfos;
}
unsigned short* rankNormalize::getRank() {
	return rank;
}

void rankNormalize::setGeneInfoByHr() {
	for (int i = 0; i < n; i++) {
		geneInfos[i].setData(hr.get_data()[i]);
		geneInfos[i].setIndices(hr.get_indices()[i]);
	}
}

bool rankNormalize::cmp1(geneInfo x, geneInfo y) {
	return x.getData() < y.getData();
}
bool rankNormalize::cmp2(geneInfo x, geneInfo y) {
	return x.getIndices() < y.getIndices();
}
void rankNormalize::sortByData(int begin,int len) {
	auto bound_cmp = bind(&rankNormalize::cmp1, this, _1, _2);
	sort(geneInfos + begin, geneInfos + begin + len, bound_cmp);
}
void rankNormalize::sortByIndices(int begin,int len) {
	auto bound_cmp = bind(&rankNormalize::cmp2, this, _1, _2);
	sort(geneInfos + begin, geneInfos + begin + len, bound_cmp);
}
void rankNormalize::ranks(int nt) {
	int m = hr.get_cell_count();
	thread *threads;
	threads = new thread[nt];
	for (int k = 0; k < m / nt + 1; k++) {
		for (int i = 0; i < min(nt, m - nt * k); i++) {
			//cout << k * nt + i << endl;
			threads[i] = thread(&rankNormalize::ranksThread, this, k*nt + i+1);
		}

		for (int i = 0; i < min(nt, m - nt * k); i++) {
			threads[i].join();
		}
	}
	delete[] threads;
}
void rankNormalize::ranksThread(int i) {
	/*int m = hr.get_cell_count();
	for (int i = 1; i <= m; i++) {*/
	//mut.lock();
		long long begin = hr.get_indptr()[i - 1];
		long long end = hr.get_indptr()[i];
		int len = end - begin;

		//cout << begin << " " << len << endl;
		//cin.get();
		//cin.get();
		sortByData(begin, len);
		
		//cout << endl;
		long long k = hr.get_gene_count() - len +1 ;
		//cout <<k<<" "<< len << endl;
		for (int j = begin; j < end; j++) {
			int cnt = 1;
			int b = j;
			while (j + 1 < end && geneInfos[j].getData() == geneInfos[j + 1].getData()) {
				cnt++;
				j++;
			}
			unsigned short value = (unsigned short)((k + (b - begin) + k + (j - begin)) / 2);
			//cout << k + b << " " << k + j << endl;
			for (int p = b; p <= j; p++) {
				geneInfos[p].setRank(value);
			}
			//cout << endl;
			
			
			/*for (int h = b; h <= j; h++) {
				cout << h << " " << geneInfos[h].getData() << " " ;
				printf("%.1f \n", geneInfos[h].getRank());
			}*/
		}
		
		sortByIndices(begin, len);
		for (int j = begin; j < end; j++) {
			rank[j] = geneInfos[j].getRank();
			//printf("%d %.1f\n", geneInfos[j].getIndices(),rank[j]);
		}
		//cout << endl;
		//mut.unlock();
	//}
}
void rankNormalize::print() {
	for (int i = 0; i < n; i++) {
		cout <<i<<" "<< rank[i] << endl;
	}
}
