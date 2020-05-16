#include "stdio.h"
#include "iostream"
#include "fstream"
#include "stdlib.h"
#include "math.h"
#include <vector>
#include <algorithm>


#define Size 256
using namespace std;

unsigned char Sample_img[Size][Size];
unsigned char Original_img[Size][Size];

class Image{
	private:
		unsigned char ImageData[Size][Size];
	public:
		Image();

        
        void low_pass(double weight){
            int filter_size = 3;
            int length = (filter_size-1)/2;
            double tmp[Size][Size];
            double mask[3][3] = {{1,weight,1}, {weight,pow(weight,2),weight}, {1,weight,1}};
                    
            for(int i=0; i<Size; i++){
                for(int j=0; j<Size; j++){
                    double convolution_val = 0;
                    for(int k=0; k<filter_size; k++){
                        if( (i-length+k) < 0 || (i-length+k) >= Size) continue;
                        for(int l=0; l<filter_size; l++){
                            if( (j-length+l) < 0 || (j-length+l) >= Size) continue;
                            convolution_val += mask[k][l] * ImageData[i-length+k][j-length+l];
                        }
                    }
                    tmp[i][j] = convolution_val;
                }
            }
            
            for(int i = 0; i < Size; i++)
				for(int j = 0; j < Size; j++)
					ImageData[i][j] = round(tmp[i][j] / pow((weight+2),2) );
                    //ImageData[i][j] = round(tmp[i][j]);

        }        

		void median_filter(int filter_sz){
			int length = (filter_sz-1)/2;
			double tmp[Size][Size];
			for(int i = 0; i < Size; i++)
				for(int j = 0; j < Size; j++){
					vector<int> vec;
					int cnt = 0;
					for(int k = i-length; k <= i+length; k++){
						if(k < 0 || k >= Size) continue;
						for(int l = j-length; l <= j+length; l++){
							if(l < 0 || l >= Size) continue;
							vec.push_back(ImageData[k][l]);
							cnt++;
						}
					}
					sort(vec.begin(), vec.end());
					tmp[i][j] = vec[round(cnt/2)];
				}
			for(int i = 0; i < Size; i++)
				for(int j = 0; j < Size; j++)
					ImageData[i][j] = tmp[i][j];
		}

		void PSNR(){
			double mse = 0;
			for(int i = 0; i < Size; i++)
				for(int j = 0; j < Size; j++)
					mse += pow(Sample_img[i][j]-ImageData[i][j], 2);
			//printf("MSE: %lf\n", mse);
			mse /= Size*Size;
			//printf("PSNR: %lfdb\n", 10*log10(255*255/mse));
            cout << "PSNR: " << 10*log10(255*255/mse) << endl;
		}
		
		void output(char out[100]){
			FILE *fout;
			if (!(fout=fopen(out,"wb"))){
				cout<<"Cannot open file!"<<endl;
				exit(1);
			}
			fwrite(ImageData, sizeof(unsigned char), Size*Size, fout);
			fclose(fout);
		}

};
Image::Image(){
	for(int i = 0; i < Size; i++)
		for(int j = 0; j < Size; j++)
			ImageData[i][j] = Original_img[i][j];
}

int main(int argc, char **argv)
{
	FILE *fin;
	if (!(fin=fopen(argv[1],"rb"))){
		cout<<"Cannot open file!"<<endl;
		exit(1);
	}
	fread(Sample_img, sizeof(unsigned char), Size*Size, fin);
	fclose(fin);
    
    
	if (!(fin=fopen(argv[2],"rb"))){
		cout<<"Cannot open file!"<<endl;
		exit(1);
	}
	fread(Original_img, sizeof(unsigned char), Size*Size, fin);
	fclose(fin);

    
    Image S4, N1;
        
    S4.PSNR();
    N1.low_pass(4);
    N1.PSNR();
    N1.output(argv[4]);
	
	if (!(fin=fopen(argv[3],"rb"))){
		cout<<"Cannot open file!"<<endl;
		exit(1);
	}
	fread(Original_img, sizeof(unsigned char), Size*Size, fin);
	fclose(fin);
	Image S5, N2;

    S5.PSNR();
	N2.median_filter(3);
    N2.PSNR();
    N2.output(argv[5]);
    
	return 0;
}