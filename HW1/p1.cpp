#include "stdio.h"
#include "iostream"
#include "fstream"
#include "stdlib.h"
#include "math.h"

#define Size 256
using namespace std;

unsigned char Original_img[Size][Size];

class Image{
	private:
		unsigned char ImageData[Size][Size];
	public:
		Image();
		void decrease_brightness(int factor){
			for(int i = 0; i < Size; i++)
				for(int j = 0; j < Size; j++)
					ImageData[i][j] /= factor;
		}

		void hist_equal(){
			int cnt[Size] = {0};
			float cdf[Size] = {0};
			for(int i = 0; i < Size; i++)
				for(int j = 0; j < Size; j++)
					cnt[ImageData[i][j]]++;
			cdf[0] = float(cnt[0]) / (Size*Size);
			for(int i = 1; i < Size; i++)
				cdf[i] = cdf[i-1] + float(cnt[i]) / (Size*Size);
			for(int i = 0; i < Size; i++)
				for(int j = 0; j < Size; j++)
					ImageData[i][j] = char(round(cdf[ImageData[i][j]]*255));
		}

		void local_hist_equal(int winsize){
			int length = (winsize-1)/2;
			unsigned char tmp[Size][Size];
			for(int i = 0; i < Size; i++){
				for(int j = 0; j < Size; j++){
					int rank = 0;
					int cnt = 0;
					for(int k = i-length; k <= i+length; k++){
						if(k < 0 || k >= Size)	continue;
						for(int l = j-length; l <= j+length; l++){
							if(l < 0 || l >= Size)	continue;
							cnt++;
							if(ImageData[i][j] > ImageData[k][l])	rank++;
						}
					}
					tmp[i][j] = round(rank*255 / cnt);
				}
			}
			for(int i = 0; i < Size; i++)
				for(int j = 0; j < Size; j++)
					ImageData[i][j] = tmp[i][j];
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

        void out_hist(char out[100]){
            //count hist
            int cnt[Size] = {0};
			for(int i = 0; i < Size; i++)
				for(int j = 0; j < Size; j++)
					cnt[ImageData[i][j]]++;
            //write txt file
            ofstream myfile (out);
            if (myfile.is_open())
            {
                myfile << "value,count" << endl;
                for(int i=0; i<Size; i++){
                    myfile << i << "," << cnt[i] << endl;
                }
                myfile.close();
            }
            else {
                cout << "Unable to open file";
                exit(1);
            }
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
	fread(Original_img, sizeof(unsigned char), Size*Size, fin);
	fclose(fin);

	Image D, E, Hd, He, Ld, Le;
    
	D.decrease_brightness(2);
	D.output(argv[2]);
	E.decrease_brightness(3);
	E.output(argv[3]);
    Hd.decrease_brightness(2);
    Hd.hist_equal();
    Hd.output(argv[4]);
    He.decrease_brightness(3);
    He.hist_equal();
    He.output(argv[5]);
    Ld.decrease_brightness(2);
    Ld.local_hist_equal(31);
    Ld.output(argv[6]);
    Le.decrease_brightness(3);
    Le.local_hist_equal(31);
    Le.output(argv[7]);

    

	return 0;
}