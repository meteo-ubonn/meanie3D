/* The MIT License (MIT)
 * 
 * (c) Jürgen Simon 2014 (juergen.simon@uni-bonn.de)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef M3D_WEIGHT_FUNCTION_FACTORY_IMPL_H
#define	M3D_WEIGHT_FUNCTION_FACTORY_IMPL_H

#include "weight_function.h"
#include "oase_weights.h"
#include "ci_weights.h"
#include "default_weights.h"
#include "inverse_default.h"
#include "exp10_weight.h"

#include "weight_function_factory.h"

#include <meanie3D/clustering/detection.h>

namespace m3D {

    template <typename T>
    WeightFunction<T> *
    WeightFunctionFactory<T>::create(const detection_params_t<T> &params,
            const detection_context_t<T>& ctx)
    {
        if (params.verbosity > VerbositySilent) {
            cout << endl << "Constructing " << params.weight_function_name << " weight function ...";
            start_timer();
        }

        WeightFunction<T> *weight_function = NULL;
        if (params.weight_function_name == "oase-ci") {
            weight_function = new OASECIWeightFunction<T>(
                    ctx.fs,
                    params.filename,
                    ctx.bandwidth,
                    params.ci_protocluster_scale,
                    params.ci_protocluster_min_size,
                    params.ci_comparison_file,
                    params.ci_comparison_protocluster_file,
                    params.ci_satellite_only,
                    params.ci_use_walker_mecikalski);
        }
        else if (params.weight_function_name == "oase") 
        {
            weight_function = new OASEWeightFunction<T>(ctx.fs, ctx.data_store, ctx.bandwidth);
        }
        else if (params.weight_function_name == "inverse") 
        {
            if (ctx.sf != NULL) {
                weight_function = new InverseDefaultWeightFunction<T>(ctx.fs,
                        ctx.data_store,
                        ctx.sf->get_filtered_min(),
                        ctx.sf->get_filtered_max());
            } else {
                weight_function = new InverseDefaultWeightFunction<T>(
                        ctx.fs,
                        ctx.data_store,
                        params.lower_thresholds,
                        params.upper_thresholds);

            }
        }
        else if (params.weight_function_name == "pow10") 
        {
            weight_function = new EXP10WeightFunction<T>(ctx.fs, ctx.data_store);
        }
        else 
        {
            // Default 
            if (ctx.sf != NULL) {
                weight_function = new DefaultWeightFunction<T>(
                        ctx.fs,
                        ctx.data_store,
                        ctx.sf->get_filtered_min(),
                        ctx.sf->get_filtered_max());
            } else {
                weight_function = new DefaultWeightFunction<T>(ctx.fs);
            }
        }

        return weight_function;
    }
}

#endif	/* M3D_WEIGHT_FUNCTION_FACTORY_IMPL_H */

