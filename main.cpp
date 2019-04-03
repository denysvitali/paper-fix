#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>

#include <iostream>
#include <queue>
#include <set>
#include <vector>
#include <algorithm>

using namespace cv;
using namespace std;

const char* wndname = "Paper Fix";
bool ROTATE_FLAG = false;

static double calculateArea(const vector<Point>& rectangle) {
  Point prevPoint;
  double area1, area2;

  double ab = sqrt(pow(rectangle[0].x - rectangle[1].x,2) + pow(rectangle[0].y - rectangle[1].y,2));
  double bc = sqrt(pow(rectangle[1].x - rectangle[2].x,2) + pow(rectangle[1].y - rectangle[2].y,2));
  area1 = (ab + bc) / 2.0;

  double cd = sqrt(pow(rectangle[2].x - rectangle[3].x,2) + pow(rectangle[2].y - rectangle[3].y, 2));
  double ad = sqrt(pow(rectangle[3].x - rectangle[0].x, 2) + pow(rectangle[3].y - rectangle[0].y,2));
  area2 = (cd + ad) / 2.0;

  return area1 + area2;
}


class Rectangle {
  public:
    Rectangle(vector<Point> v) : m_v{v} {
      /*for(Point& p: v){
        cout << p << ",";
      }
      cout << endl;*/
      computeArea();
    }
    double area() const {
      return m_area;
    }
    const Point& operator [] (const int index){
      return m_v[index];
    }

    size_t size() const {
      return m_v.size();
    }

    vector<Point> vec() const {
      return m_v;
    }

    bool equals(const Rectangle& o) const {
      for(size_t i=0; i<o.m_v.size(); i++){
        if(m_v[i].x != o.m_v[i].x){
          return false;
        }

        if(m_v[i].y != o.m_v[i].y){
          return false;
        }
      }
      cout << "Equals! " << endl;
      return true;
    }
  private:
    vector<Point> m_v;
    double m_area = 0.0;

    void computeArea(){
      m_area = calculateArea(m_v);
    }
};

void usage(int argc, char* argv[]){
  cout << "Usage: ";
  cout << argv[0] << " [options] input.jpg output.jpg" << endl;
  cout << "Options: " << endl;
  cout << "-r / --rotate\t\tRotates the final image by 180° (upside down)" << endl;
}

static double angle( Point pt1, Point pt2, Point pt0 )
{
    double dx1 = pt1.x - pt0.x;
    double dy1 = pt1.y - pt0.y;
    double dx2 = pt2.x - pt0.x;
    double dy2 = pt2.y - pt0.y;
    return (dx1*dx2 + dy1*dy2)/sqrt((dx1*dx1 + dy1*dy1)*(dx2*dx2 + dy2*dy2) + 1e-10);
}

void find_rectangles(Mat& image, vector<vector<Point> >& rectangles)
{
    // blur will enhance edge detection
    Mat blurred(image);
    medianBlur(image, blurred, 9);

    Mat gray0(blurred.size(), CV_8U), gray;
    vector<vector<Point>> contours;

    // find rectangles in every color plane of the image
    for (int c = 0; c < 3; c++)
    {
        int ch[] = {c, 0};
        mixChannels(&blurred, 1, &gray0, 1, ch, 1);

        // try several threshold levels
        const int threshold_level = 2;
        for (int l = 0; l < threshold_level; l++)
        {
            // Use Canny instead of zero threshold level!
            // Canny helps to catch rectangles with gradient shading
            if (l == 0)
            {
                Canny(gray0, gray, 10, 20, 3); //
                // Dilate helps to remove potential holes between edge segments
                dilate(gray, gray, Mat(), Point(-1,-1));
            }
            else
            {
                    gray = gray0 >= (l+1) * 255 / threshold_level;
            }

            // Find contours and store them in a list
            findContours(gray, contours, RETR_LIST, CHAIN_APPROX_SIMPLE);

            // Test contours
            vector<Point> approx;
            for (size_t i = 0; i < contours.size(); i++)
            {
                    // approximate contour with accuracy proportional
                    // to the contour perimeter
                    approxPolyDP(Mat(contours[i]), approx, arcLength(Mat(contours[i]), true)*0.02, true);

                    // Note: absolute value of an area is used because
                    // area may be positive or negative - in accordance with the
                    // contour orientation
                    if (approx.size() == 4 &&
                            fabs(contourArea(Mat(approx))) > 1000 &&
                            isContourConvex(Mat(approx)))
                    {
                            double maxCosine = 0;

                            for (int j = 2; j < 5; j++)
                            {
                                    double cosine = fabs(angle(approx[j%4], approx[j-2], approx[j-1]));
                                    maxCosine = MAX(maxCosine, cosine);
                            }

                            if (maxCosine < 0.3)
                              rectangles.push_back(approx);
                    }
            }
        }
    }
}

/*
static int getPointIndex(vector<Point>& vp, Point2f& p){
  for(size_t i=0; i<vp.size(); i++){
    if(vp[i].x == p.x && vp[i].y == p.y){
      return i;
    }
  }
  return -1;
}
*/

template <class T> static Mat* drawRectangles(Mat& image, Mat& originalImage, set<Rectangle, T>& rectangles, int rows, int cols )
{
    auto cmp = [](Rectangle& f, Rectangle& s) { return f.area() <= s.area(); };
    priority_queue<Rectangle, std::vector<Rectangle>, decltype(cmp)> pq(cmp);
    double rotation = 0;
    int i = 0;

    cout << "Found " << rectangles.size() << " rectangles" << endl;

    for(const Rectangle& s : rectangles)
    {
        pq.push(s);
        //double area = s.area();
        //cout << "Rectangle " << i << ": Area = " << area << endl;
        i++;
    }

    int count = 0;

    vector<Scalar> colorArr = {Scalar(0,255,0),
      Scalar(0,0,255), Scalar(255,0,0)};

    while(!pq.empty() && count < 1){
      const Rectangle& s = pq.top();
      pq.pop();
      vector<Point> points = s.vec();

      const Point* p = &(points[0]);
      int n = (int) s.size();

      bool insideOrigImage = true;
      for(Point& p: points){
        if(!p.inside(Rect(0,0, cols * 1.05, rows * 1.05))){
          insideOrigImage = false;
          break;
        }
      }

      if (p-> x > 3 && p->y > 3 && insideOrigImage){
        if(rotation == 0) {
          rotation = atan2(p[0].y, p[1].y);
          //cout << "Rotation: " << rotation << endl;
        }
        polylines(image, &p, &n, 1, true, colorArr[count%3], 3, LINE_AA);

        Rect boundRect = boundingRect(points);
        vector<Point2f> quad_pts, square_pts;

        vector<Point> points_clone = points;
        // cout << "Points: " << points_clone << endl;

        Point2f topLeft, topRight, bottomLeft, bottomRight;
        sort(points_clone.begin(), points_clone.end(), [](const Point& lhs, const Point& rhs){
          return lhs.y < rhs.y;
        });


        if(points_clone[0].x > points_clone[1].x){
           topLeft = points_clone[1];
           topRight = points_clone[0];
        } else {
           topLeft = points_clone[0];
           topRight = points_clone[1];
        }

        if(points_clone[2].x > points_clone[3].x){
           bottomLeft = points_clone[3];
           bottomRight = points_clone[2];
        } else {
           bottomLeft = points_clone[2];
           bottomRight = points_clone[3];
        }
        
        

        quad_pts.push_back(topLeft);
        quad_pts.push_back(bottomLeft);
        quad_pts.push_back(topRight);
        quad_pts.push_back(bottomRight);

        square_pts.push_back(Point2f(boundRect.x, boundRect.y));
        square_pts.push_back(Point2f(boundRect.x, boundRect.y + boundRect.height));
        square_pts.push_back(Point2f(boundRect.x + boundRect.width, boundRect.y));
        square_pts.push_back(Point2f(boundRect.x + boundRect.width, boundRect.y + boundRect.height));

        rectangle(image,boundRect,Scalar(255,0,0),1,8,0);

        Mat transmtx = getPerspectiveTransform(quad_pts,square_pts);
        Mat transformed = Mat::zeros(image.rows, image.cols, CV_8UC3);
        warpPerspective(originalImage, transformed, transmtx, image.size());
        //imwrite("transformed.jpg", transformed);

        Mat* ROI = new Mat(transformed, Rect(boundRect.x, boundRect.y, boundRect.width, boundRect.height));
        //imwrite("cropped.jpg", *ROI);

        return ROI;
      }
    }
    return nullptr;
}

static Mat* expandImage( Mat& src ){
  int top = (int) (0.05*src.rows); 
  int bottom = (int) (0.05*src.rows);
  int left = (int) (0.05*src.cols); 
  int right = (int) (0.05*src.cols);

  Mat* dst = new Mat{};
  copyMakeBorder(src, *dst, top, bottom, left, right, BORDER_CONSTANT, 0);

  return dst;
}

double distance(Point& p1, Point& p2){
  return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

void rotate(cv::Mat& src, double angle, cv::Mat& dst){
    cv::Point2f ptCp(src.cols*0.5, src.rows*0.5);
    cv::Mat M = cv::getRotationMatrix2D(ptCp, angle, 1.0);
    cv::warpAffine(src, dst, M, src.size(), cv::INTER_CUBIC); //Nearest is too rough, 
}

int main(int argc, char* argv[]){

  int options_pos = 1;
  for(int i=1; i<argc; i++){
    if(strcmp(argv[i], "-r") == 0 || strcmp(argv[i], "--rotate") == 0){
      ROTATE_FLAG=true;
      options_pos++;
    }
  }

  Mat image;
  if(options_pos + 2 > argc){
    cerr << "Invalid input" << endl;
    usage(argc, argv);
    return 1;
  }

  image = imread(argv[options_pos]);

  if(!image.data){
    cerr << "No image data." << endl;
    return 1;
  }

  //Mat gray_image;
  //cvtColor( image, gray_image, CV_BGR2GRAY );
  
  vector<vector<Point>> rectangles;
  Mat expandedImage = *expandImage(image);

  Mat originalImageExpanded;
  expandedImage.copyTo(originalImageExpanded);

  find_rectangles(expandedImage, rectangles);

  auto cmp = [](Rectangle a, Rectangle b) { return !a.equals(b); };
  set<Rectangle, decltype(cmp)> theRectangles(cmp);

  Point topLeft {(int) (0.05*image.cols), (int) (0.05*image.rows)};
  Point bottomLeft {(int) (0.05*image.cols), (int) (1.05*image.rows)};

  //cout << "Top: " << topLeft << ", BottomLeft: " << bottomLeft << endl;

  for(vector<Point> p : rectangles){
    //cout << p << endl;
    if(distance(p[0],topLeft) < 5 && distance(p[1],bottomLeft) < 5){
      cout << "Ignoring this rectangle because it's basically a rectangle containing our picture." << endl;
      continue;
    }
    theRectangles.insert(Rectangle(p));
  }
  Mat* finalImage = drawRectangles(expandedImage, originalImageExpanded, theRectangles, image.rows, image.cols);

  if(ROTATE_FLAG){
    // Rotate by 180°
    Mat rotatedFinal;
    rotate(*finalImage, rotatedFinal, 1);
    imwrite(argv[options_pos + 1], rotatedFinal);
  } else {
    Mat outputImage = *finalImage;
    imwrite(argv[options_pos + 1], outputImage);
  }

  delete finalImage;
}
