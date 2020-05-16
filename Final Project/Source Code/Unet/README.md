### Steps for running Unet 
1. put these 3 python files (main.py, data.py, model.py) under the same folder
2. create new folder and name it "data"
3. put the training images into "data/train/image" folder
4. manually annotate images and put them into "data/train/label" folder
5. put testing images into "data/test" folder
6. execute main.py
7. the results will be saved in "data/test" folder and named as "_predict.png"

### p.s. 因為我們所用的影像資料(source images)是來自醫生，不能隨意外流，所以沒有附上