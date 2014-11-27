#ifndef M3D_OPENCV_UTILS_H
#define M3D_OPENCV_UTILS_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#if WITH_OPENCV
#include <opencv2/opencv.hpp>

namespace m3D { namespace utils { 

    #define VISUALIZE_FLOW_FIELD 1

    template <typename T>
    cv::Mat matrix_from_multiaray(MultiArray<T> *array,
                                  T min,
                                  T max,
                                  bool invert = false)
    {

        vector<size_t> dimension_sizes = array->get_dimensions();

        if (dimension_sizes.size() != 2)
        {
            cerr << "FATAL:this method is only supported for 2D data" << endl;
            exit(EXIT_FAILURE);
        }

        using namespace cv;

        // transform the data into a grey image with values
        // between 0 (valid_min) and 255 (valid_max)

        Mat image;
        image.create( dimension_sizes[0], dimension_sizes[1], CV_8U);

        for (size_t i=0; i<dimension_sizes[0]; i++)
        {
            for (size_t j=0; j<dimension_sizes[1]; j++)
            {
                vector<int> gridpoint(2);
                gridpoint[1] = i;
                gridpoint[0] = j;

                bool isValid = false;

                T value = array->get(gridpoint);

                T scaled_value = (256.0/(max-min)) * (value - min);

                unsigned char grayValue = (unsigned char) round(scaled_value);

                if (invert)
                {
                    grayValue = 255.0 - grayValue;
                }

                image.at<uchar>( dimension_sizes[1] - j - 1, i, 0 ) = grayValue;
            }
        }

        return image;
    }

    template <typename T>
    cv::Mat matrix_from_variable(DataStore<T> *dataStore,
                                 size_t variableIndex,
                                 T lower_treshold = std::numeric_limits<T>::min(),
                                 T upper_treshold = std::numeric_limits<T>::max(),
                                 bool invert = false,
                                 int timeIndex=-1)
    {
        vector<size_t> dimension_sizes = dataStore->get_dimension_sizes();

        if (dimension_sizes.size() != 2)
        {
            cerr << "FATAL:this method is only supported for 2D data" << endl;
            exit(EXIT_FAILURE);
        }

        using namespace cv;

        // transform the data into a grey image with values
        // between 0 (valid_min) and 255 (valid_max)

        T min = dataStore->min(variableIndex);
        T max = dataStore->max(variableIndex);

        Mat image;
        image.create( dimension_sizes[0], dimension_sizes[1], CV_8U);

        for (size_t i=0; i<dimension_sizes[0]; i++)
        {
            for (size_t j=0; j<dimension_sizes[1]; j++)
            {
                vector<int> gridpoint(2);
                gridpoint[1] = i;
                gridpoint[0] = j;

                bool isValid = false;

                T value = dataStore->get(variableIndex,gridpoint,isValid);

                if (isValid && value <= upper_treshold && value >= lower_treshold)
                {
                    T scaled_value = (256.0/(max-min)) * (value - min);

                    unsigned char grayValue = (unsigned char) round(scaled_value);

                    if (invert)
                    {
                        grayValue = 255.0 - grayValue;
                    }

                    image.at<uchar>( dimension_sizes[1] - j - 1, i, 0 ) = grayValue;
                }
            }
        }

        return image;
    }

    /** Helper method. Visualizes a flow field
     */
    static inline void
    drawOptFlowMap(const cv::Mat& flow, cv::Mat& cflowmap, int step, double, const cv::Scalar& color)
    {
        using namespace cv;

        for(int y = 0; y < cflowmap.rows; y += step)
            for(int x = 0; x < cflowmap.cols; x += step)
            {
                const Point2f& fxy = flow.at<Point2f>(y, x);
                line(cflowmap, cv::Point(x,y), cv::Point(cvRound(x+fxy.x), cvRound(y+fxy.y)),
                     color);
                circle(cflowmap, cv::Point(x,y), 2, color, -1);
            }
    }

    /** Evaluate optical flow from the two data stores in the given variable
     * and shift the previous data along the flow vectors to the current 
     * position. 
     *
     * Flow estimation is done with OpenCV, using the algorithm implementing
     * the method described in this paper:
     *
     * Farnebäck, Gunnar. “Two-Frame Motion Estimation Based on Polynomial Expansion.” 
     * In Image Analysis, 363–70. Springer, 2003. 
     * http://link.springer.com/chapter/10.1007/3-540-45103-X_50. 
     *
     * This method has the advantage of providing a dense vector field, which 
     * makes morphing easier. 
     * 
     * @param filename of the current netcdf file
     * @param filename of the previous netcdf file
     * @param coordinate system
     * @param variable to use for shifting
     * @param lower threshold for that variable
     * @param upper threshold for that variable
     * @param time index
     * @param if this flag is set, the values are inverted for the flow calcualation.
     *
     * @return an instance of NetCDFDataStore where all variables have been
     * shifted along the flow field calculated from the specified variable
     */
    template <typename T>
    NetCDFDataStore<T> *
    shifted_store_from_flow_of_variable(const std::string &current,
                                        const std::string &previous,
                                        const CoordinateSystem<T> *coord_system,
                                        const std::vector<std::string> variables,
                                        const size_t &flow_variable_index,
                                        const T &lower_threshold = std::numeric_limits<T>::min(),
                                        const T &upper_threshold = std::numeric_limits<T>::max(),
                                        const int time_index = -1,
                                        bool invert_values = false)
    {
        using namespace cv;

        if (coord_system->rank() != 2)
        {
            cerr << "FATAL:'shifted_store_from_flow_of_variable' is only supported for 2D" << endl;
            exit(EXIT_FAILURE);
        }

        Mat previous_image_raw, previous_image, current_image_raw, current_image, flow, cflow;

        // Read from NetCDF
        NetCDFDataStore<T> *curr_data = new NetCDFDataStore<T>(current,coord_system,variables,time_index);
        current_image_raw = matrix_from_variable(curr_data, flow_variable_index, lower_threshold, upper_threshold, false);

        //increase the contrast (5x)
        current_image_raw.convertTo(current_image, -1, 2, 0);
        current_image_raw.release();

        // Read from NetCDF

        NetCDFDataStore<T> *prev_data = new NetCDFDataStore<T>(previous,coord_system,variables,time_index);
        previous_image_raw = matrix_from_variable(prev_data, flow_variable_index, lower_threshold, upper_threshold, false);

        //increase the contrast (5x)
        previous_image_raw.convertTo(previous_image, -1, 2, 0);
        previous_image_raw.release();

        cout << "Calculating optical flow ...";
        start_timer();
        calcOpticalFlowFarneback(previous_image, current_image, flow, 0.5, 3, 20/*km*/, 3, 10, 1.2, 0);
//         calcOpticalFlowSF(previous_image, current_image, flow, 3, 15, 15);
//        calcOpticalFlowSF(previous_image, current_image, flow, 3, 2, 4, 4.1, 25.5, 18, 55.0, 25.5, 0.35, 18, 55.0, 25.5, 10);

//        Mat transform = estimateRigidTransform(previous_image,current_image,true);

        cout << " done ("<<stop_timer()<<"s)." << endl;


#if VISUALIZE_FLOW_FIELD

        std::string prev_window_name = variables[flow_variable_index] + " (previous)";
        std::string curr_window_name = variables[flow_variable_index] + " (current)";
        std::string flow_window_name = variables[flow_variable_index] + " (flow)";

        namedWindow(curr_window_name, 1);
        namedWindow(prev_window_name, 1);
        namedWindow(flow_window_name, 1);

        imshow(prev_window_name, previous_image);
        imshow(curr_window_name, current_image);

        cvtColor(previous_image, cflow, COLOR_GRAY2BGR);
        drawOptFlowMap(flow, cflow, 10, 1.5, Scalar(0, 255, 0));
        imshow(flow_window_name, cflow);
        cflow.release();

        cvWaitKey(0);

        cvDestroyWindow(prev_window_name.c_str());
        cvDestroyWindow(curr_window_name.c_str());

        previous_image.release();
        current_image.release();
#endif

        // Iterate over the data store's variables and
        // morph them accordingly

        for (size_t vi=0; vi < variables.size(); vi++)
        {
            T NOT_SET = prev_data->min(vi);

            vector<size_t> dims = coord_system->get_dimension_sizes();

            // Initialize new array with a value for 'not found'

            MultiArray<T> *dest = new MultiArrayBlitz<T>(dims,NOT_SET);

            for(int y = 0; y < flow.rows; y += 1)
            {
                for(int x = 0; x < flow.cols; x += 1)
                {
                    vector<int> source;
                    source.push_back(y);
                    source.push_back(x);

                    bool is_valid = false;
                    T value = prev_data->get(vi,source,is_valid);

                    if (is_valid && value >= lower_threshold   && value <= upper_threshold)
                    {
                        const Point2f& fxy = flow.at<Point2f>(y, x);
                        int xdest = x + fxy.x;
                        int ydest = y + fxy.y;

                        if (xdest >= 0 && xdest < dims[1] && ydest >= 0 && ydest < dims[0])
                        {
                            // copy source to dest, provided the destination
                            // is inside the bounds

                            vector<int> gp;
                            gp.push_back(ydest);
                            gp.push_back(xdest);

                            dest->set(gp,value);
                        }
                    }
                }
            }

            // TODO: make a second pass and interpolate missing values
            // by averaging neighbours

#if VISUALIZE_FLOW_FIELD
            namedWindow(variables[vi], 1);
            Mat before = matrix_from_variable(prev_data, vi);
            imshow(variables[vi], before);
            cvWaitKey(0);
            before.release();
#endif
            prev_data->set_data(vi,dest);

#if VISUALIZE_FLOW_FIELD
            Mat after = matrix_from_variable(prev_data, vi);
            imshow(variables[vi], after);
            cvWaitKey(0);
            after.release();
            cvDestroyWindow(variables[vi].c_str());
#endif
        }

        flow.release();

#if VISUALIZE_FLOW_FIELD
        cvDestroyWindow(flow_window_name.c_str());
#endif

        delete curr_data;

        return prev_data;
    }

    template<typename T>
    void
    display_variable(NetCDFDataStore<T> *ds, int var_index)
    {
        using namespace cv;

        Mat image = matrix_from_variable(ds, var_index);

        std::string window_name = ds->variable_names()[var_index];

        namedWindow(window_name, 1);

        imshow(window_name, image);

        cvWaitKey(0);

        cvDestroyWindow(window_name.c_str());

        image.release();
    }

    template<typename T>
    void
    display_array(MultiArray<T> *ds, T min, T max, bool invert=false)
    {
        using namespace cv;

        Mat image = matrix_from_multiaray(ds,min,max,invert);

        std::string window_name = "raw array output";

        namedWindow(window_name, 1);

        imshow(window_name, image);

        cvWaitKey(0);

        cvDestroyWindow(window_name.c_str());

        image.release();
    }
}}
#endif
#endif
