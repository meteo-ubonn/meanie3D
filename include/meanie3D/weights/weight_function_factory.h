/* 
 * File:   weight_function_factory.h
 * Author: simon
 *
 * Created on January 18, 2016, 9:07 PM
 */

#ifndef WEIGHT_FUNCTION_FACTORY_H
#define	WEIGHT_FUNCTION_FACTORY_H

#include <meanie3D/weights/weight_function.h>
#include <meanie3D/clustering/detection.h>

namespace m3D {
    
    /** 
     * 
     */
    template <typename T>
    class WeightFunctionFactory {

        public:
            
            /**
             * Create a weight function from information in a detection
             * parameter and context set. 
             * 
             * @param params
             * @param ctx
             * @return 
             */
            static 
            WeightFunction<T> *create(const detection_params_t<T> &params, 
                const detection_context_t<T>& ctx);
    };
    
}

#endif	/* WEIGHT_FUNCTION_FACTORY_H */

