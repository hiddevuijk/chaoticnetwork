#ifndef GUARD_generate_matrix_h
#define GUARD_generate_matrix_h


std::vector<std::vector<int> > connectivity_matrix(int N1, int N2, int K, Ranq1& r)
{
	std::vector<int> row_vec(N2,0);
	std::vector<std::vector<int> > C(N1,row_vec);

	for (int i=0;i<N1;i++) {
		std::vector<int> row(N2,0);
		for(int n=0;n<K;n++) row[n] = 1;
		for(int n=0;n<N2;n++) {
			int r_int = r.int32() % (N2-1);
			int a = row[n];
			row[n] = row[r_int];
			row[r_int] = a;
		}
		for(int j=0;j<N2;j++) C[i][j] = row[j]; 
	}

	return C;
}

int dotproduct(const std::vector<int>& v1, const std::vector<int>& v2,const int& N)
{
	int sum = 0;
	for(int i=0;i<N;i++)
		sum += v1[i]*v2[i];
	return sum;
}

#endif

