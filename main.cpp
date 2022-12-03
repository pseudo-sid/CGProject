#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
using namespace std;
using namespace cv;

vector<pair<int, int>> coordinates;
void CallBackFunc(int event, int x, int y, int flags, void* userdata)
{
    if  ( event == EVENT_LBUTTONDOWN )
    {
        cout << "Left button of the mouse is clicked - position (" << x << ", " << y << ")" << endl;
        pair<int, int> back = coordinates.back();
        if(back.first != x or back.second != y)
            coordinates.push_back({x, y});
    }
    else if  ( event == EVENT_RBUTTONDOWN )
    {
        cout << "Right button of the mouse is clicked - position (" << x << ", " << y << ")" << endl;
    }
    else if  ( event == EVENT_MBUTTONDOWN )
    {
        cout << "Middle button of the mouse is clicked - position (" << x << ", " << y << ")" << endl;
    }
//    else if ( event == EVENT_MOUSEMOVE )
//    {
//        cout << "Mouse move over the window - position (" << x << ", " << y << ")" << endl;
//
//    }
}


int main() {
    Mat img = imread("tree.jpeg");
    if(img.empty()){
        cout << "Image couldn't be loaded\n";
        return -1;
    }

    namedWindow("Annotate crown", 1);
    setMouseCallback("Annotate crown", CallBackFunc, NULL);

    imshow("Annotate crown", img);
    waitKey(0);

    vector<pair<int, int>> crown_coords(coordinates.begin(), coordinates.end());

    coordinates.clear();
    namedWindow("Annotate bark", 1);
    setMouseCallback("Annotate bark", CallBackFunc, NULL);

    imshow("Annotate bark", img);
    waitKey(0);

    vector<pair<int, int>> bark_coords(coordinates.begin(), coordinates.end());
    coordinates.clear();

    cout << crown_coords.size() << ", " << bark_coords.size() <<"\n";

    return 0;
}
