
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <complex>
#include <fftw3.h>
#include <math.h>

#define n 1000000
#define m 150000

#define degree_1 (m + n -2)
#define degree_2 (n - 1)
#define N (degree_1 + degree_2 + 1)
	

int main(int argc,char* argv[])
{  
	//define every parameter

        FILE *f_random;
        FILE *f_output_data_prod;

        double *in_1,*in_2,*in;
        fftw_complex *out_1,*out_2,*out;
        fftw_plan p_1, p_2, p;        

        unsigned char* input_data_in_1=(unsigned char*)malloc((n + m - 1)*sizeof(unsigned char));
        memset(input_data_in_1,0,(n + m - 1)*sizeof(unsigned char));

        unsigned char* input_data_in_2=(unsigned char*)malloc(n*sizeof(unsigned char));
        memset(input_data_in_2,0,n*sizeof(unsigned char));

        f_random=fopen("data.bin","rb");
        fread(input_data_in_1,(n + m - 1),sizeof(unsigned char),f_random);
        fread(input_data_in_2,n,sizeof(unsigned char),f_random);
        fclose(f_random);

        int* prod=(int*)malloc(N*sizeof(int));
        memset(prod,0,N*sizeof(int));

        int* dest=(int*)malloc(m*sizeof(int));
        memset(dest,0,m*sizeof(int));

	in = (double*)fftw_malloc(N*sizeof(double));  
	in_1 = (double*)fftw_malloc(N*sizeof(double));
	in_2 = (double*)fftw_malloc(N*sizeof(double));
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	out_1 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	out_2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);

	//
	memset(in, 0, N*sizeof(double));
	memset(in_1, 0, N*sizeof(double));
	memset(in_2, 0, N*sizeof(double));

	//����Ҷ�任��׼����pΪ����Ҷ���任��׼����p_1Ϊ����Ҷ�任��׼����P_2Ϊ����Ҷ�任��׼��
	p = fftw_plan_dft_c2r_1d(N, out, in, FFTW_ESTIMATE);
	p_1 = fftw_plan_dft_r2c_1d(N, in_1, out_1, FFTW_ESTIMATE);
	p_2 = fftw_plan_dft_r2c_1d(N, in_2, out_2, FFTW_ESTIMATE);


	//��ʼ������
	//initialize in_1,in_2;
	for (int i = 0; i <= degree_1; i++)
        {
           if(input_data_in_1[i]%2 == 0)
           {
              in_1[i]=0.0;
           }
           else
           {
              in_1[i]=1.0;
           }
	}
	for (int i = 0; i < degree_2; i++)
        {
	   if(input_data_in_2[i]%2 == 0)
           {
              in_2[i]=0.0;
           }
           else
           {
              in_2[i]=1.0;
           }
	}

	//ִ�и���Ҷ�任
	fftw_execute(p_1); /* repeat as needed */
	fftw_execute(p_2); /* repeat as needed */

	
	//������Ҷ�任��Ľ���������
	//cal out
	for (int i = 0; i <= N / 2; i++){
		out[i][0] = out_1[i][0] * out_2[i][0] - out_1[i][1] * out_2[i][1];
		out[i][1] = out_1[i][0] * out_2[i][1] + out_1[i][1] * out_2[i][0];
	}

	//ִ�з��任
	fftw_execute(p);

	//������
	for (int i = 0; i<N; i++)
        {
            prod[i] = in[i] / N + 0.5;
		//	std::cerr<<in[i]<<" "<<N<<" "<<in[i] / N<<" "<<prod[i]<<std::endl;
	}

        for (int i = 0; i<m; i++)
        {
            dest[i] = prod[n -1 + i] % 2;
		//	std::cerr<<in[i]<<" "<<N<<" "<<in[i] / N<<" "<<prod[i]<<std::endl;
	}

        f_output_data_prod=fopen("result.bin","wb");
        fwrite(dest,m,sizeof(int),f_output_data_prod);
        fclose(f_output_data_prod);
        
        free(prod);
        free(dest);
   
	//�ͷ�
	//delete
	fftw_destroy_plan(p);
	fftw_destroy_plan(p_1);
	fftw_destroy_plan(p_2);
        free(input_data_in_1);
        free(input_data_in_2);
  
	return 0;

}
