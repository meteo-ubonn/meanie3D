#ifndef _M3D_TRACKING_CLASS_Impl_H_
#define _M3D_TRACKING_CLASS_Impl_H_

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <vector>
#include <netcdf>

#include <numericalrecipes/nr.h>

#include <meanie3D/types/cluster.h>
#include <meanie3D/tracking/tracking.h>
#include <meanie3D/utils/verbosity.h>

namespace m3D {
    
    using namespace utils;
    using namespace netCDF;

    template <typename T>
    typename Tracking<T>::matrix_t
    Tracking<T>::create_matrix(size_t width, size_t height)
    {
        matrix_t matrix;
        
        matrix.resize(width);
        
        for (int i=0; i<width; ++i)
        {
            matrix[i].resize(height);
        }
        
        return matrix;
    }

    
    template <typename T>
    void
    Tracking<T>::track(const ClusterList<T> &last_clusters,
                       const ClusterList<T> &current_clusters,
                       const vector<NcVar> &feature_variables,
                       const NcVar &track_variable,
                       Verbosity verbosity )
    {
        // sanity check
        
        if ( current_clusters.clusters.size() == 0 )
        {
            if ( verbosity >= VerbosityNormal )
            {
                cout << "No current clusters. Skipping tracking." << endl;
                
                return;
            }
        }
        
        
        // figure out tracking variable index
        
        int index = index_of_first( feature_variables, track_variable );
        
        assert( index >= 0 );
        
        size_t tracking_var_index = (size_t) index;
        
        
        // Valid min/max
        
        T valid_min, valid_max;
        
        track_variable.getAtt("valid_min").getValues( &valid_min );

        track_variable.getAtt("valid_max").getValues( &valid_max );
        
        
        // Find the highest id 
        
        size_t current_id = Cluster<T>::NO_ID;

        // Find out the highest id from the previous list
        
        typename Cluster<T>::list::iterator ci;
        
        typename Cluster<T>::list previous = last_clusters.clusters;
        
        typename Cluster<T>::list current = current_clusters.clusters;
        
        for ( ci = previous.begin(); ci != previous.end(); ci++ )
        {
            typename Cluster<T>::ptr c = *ci;
            
            if ( c->id > current_id )
            {
                current_id = c->id;
            }
        }
        
        current_id++;
        
        if ( verbosity >= VerbosityDetails )
        {
            cout << "Tracking: next available id is " << current_id << endl;
        }
        
        // TODO: check time difference and determine displacement restraints
        
//        Blob* oldBlob = (Blob*)[blobs objectAtIndex:0];
//        Blob* newBlob = (Blob*)[newBlobs objectAtIndex:0];
//        int deltaT = [[newBlob timestamp] timeIntervalSinceDate:[oldBlob timestamp]];
//        if (doTracking) {
//            if (deltaT<0) {
//                SPLog(@"Can't go back in time when tracking. Ignored.");
//                return;
//            }
//            else if (deltaT==0) {
//                SPLog(@"Repeated Image while tracking. Nothing to do");
//                return;
//            }
//            else if ([trackStore constrainTimesliceDifference] && deltaT>[trackStore maxTimesliceDifference]) {
//                SPLog(@"Critical time difference of %4.1d exceeded while tracking. Starting new Series",
//                      [trackStore maxTimesliceDifference]);
//                [blobs removeAllObjects];
//                [self addSet:theNewBlobs toArray:blobs tagObjects:YES];
//                [trackStore addBlobArray:blobs];
//                meanVelocity=[trackStore maxVelocity];
//                return;
//            }
//        }
//        float maxDisplacement = [trackStore maxVelocity]*deltaT;
//        SPLog(@"max velocity  constraint at %4.1f m/s at deltaT %4.0fs -> dR_max = %7.1fm",
//              [trackStore maxVelocity],deltaT,maxDisplacement);
//        
//        float maxMeanVelocityDisplacement=0;
//        if (useMeanVelocityConstraint) {
//            maxMeanVelocityDisplacement=deltaT*meanVelocity*meanVelocitySecurityPercentage;
//            SPLog(@"mean velocity constraint at %4.2f m/s, dR_max = %7.1f",
//                  meanVelocity*meanVelocitySecurityPercentage,maxMeanVelocityDisplacement);
//        }
//        SPLog(@"\n");

        
        // Get the counts
        size_t old_count = previous.size();
        size_t new_count = current.size();
        
        // Create result matrices
        matrix_t rank_correlation = create_matrix(new_count,old_count);
        matrix_t midDisplacement = create_matrix(new_count,old_count);
        matrix_t histDiff = create_matrix(new_count,old_count);
        matrix_t sum_prob = create_matrix(new_count,old_count);
        matrix_t coverOldByNew = create_matrix(new_count,old_count);
        matrix_t coverNewByOld = create_matrix(new_count,old_count);
        
        // compute correlation matrix
        
        int n,m;
        
        // check out the values for midDisplacement, histSizeDifference and kendall's tau
        
        for ( n=0; n < new_count; n++ )
        {
            typename Cluster<T>::ptr newCluster = current[n];
            
            typename Histogram<T>::ptr newHistogram = newCluster->histogram(tracking_var_index,valid_min,valid_max);
            
            for ( m=0; m < old_count; m++ )
            {
                typename Cluster<T>::ptr oldCluster = previous[m];
                
                typename Histogram<T>::ptr oldHistogram = oldCluster->histogram(tracking_var_index,valid_min,valid_max);
                
                // calculate spearman's tau
                
                rank_correlation[n][m] = newHistogram->correlate_kendall( oldHistogram );
                
                rank_correlation[n][m] = newHistogram->correlate_spearman( oldHistogram );
//                
//                
//
//                Blob* oldBlob = (Blob*)[blobs objectAtIndex:m];
//                [self correlate:oldBlob byKendallWith:newBlob tau:&tau[n][m] prob:&probrs[n][m]];
//                
//                float div = GSL_MAX((float)[newBlob histogramSize],(float)[oldBlob histogramSize]);
//                if (div==0.0) div=1.0;
//                histDiff[n][m]=abs((float)[newBlob histogramSize]-(float)[oldBlob histogramSize])/div;
//                
//                float dx = [oldBlob area]->scaleX * ([oldBlob midX]-[newBlob midX]);
//                float dy = [oldBlob area]->scaleY * ([oldBlob midY]-[newBlob midY]);
//                midDisplacement[n][m]=sqrt(dx*dx + dy*dy);
//                
//                coverNewByOld[n][m]=[self coverageOf:newBlob by:oldBlob];
//                coverOldByNew[n][m]=[self coverageOf:oldBlob by:newBlob];
            }
        }
        
        if ( verbosity >= VerbosityDetails )
        {
            cout << "Correlation table " << endl;
        
            cout << "      ";
            
            cout << fixed << setprecision(4);
            
            for ( m=0; m < old_count; m++ )
            {
                cout << m;
            }
            
            cout << endl;

            for ( n=0; n < new_count; n++ )
            {
                cout << n;
                
                for ( m=0; m < old_count; m++ )
                {
                    cout << rank_correlation[n][m];
                }
                
                cout << endl;
            }
        }

        
        // find maximum values for midD and histSD
//        int maxHistD = 0;
//        float maxMidD = 0;
//        for (n=0;n<[newBlobs count];n++)
//        {
//            for (m=0;m<[blobs count];m++)
//            {
//                if (histDiff[n][m]>maxHistD)
//                    maxHistD = histDiff[n][m];
//                
//                if (midDisplacement[n][m]>maxMidD)
//                    maxMidD = midDisplacement[n][m];
//            }
//        }
//        if (maxMidD==0) maxMidD=1;
//        if (maxHistD==0) maxHistD=1;
//        
//        // normalise with maximum possible distance
//        // maxMidD = sqrt( 2*200*200 );
//        
//        // calculate the final probability
//        SPLog(@"-- Correlation Table --\n");
//        for (n=0;n<[newBlobs count];n++)
//        {
//            Blob* newBlob = [newBlobs objectAtIndex:n];
//            SPLog(@"\nCorrelating new Blob #%d (histogram size %d, [%3d,%3d])",
//                  n,(int)[newBlob histogramSize],[newBlob midX],[newBlob midY]);
//            
//            for (m=0;m<[blobs count];m++)
//            {
//                Blob* cmpBlob = (Blob*)[blobs objectAtIndex:m];
//                float prob_r = WEIGHT_R * gsl_sf_erfc(midDisplacement[n][m]/maxMidD);
//                float prob_h = WEIGHT_H * gsl_sf_erfc(histDiff[n][m]/maxHistD);
//                float prob_t = WEIGHT_T * tau[n][m];
//                sum_prob[n][m] = prob_t + prob_r + prob_h;
//                
//                SPLog(@"\t<ID#%ld>:\tdeltaR=%4.1f (%f)\tdH=%5.4f (%f)\ttau=%5.4f\tsum=%5.4f\tcovON=%3.2f\tcovNO=%3.2f",
//                      [cmpBlob tag],
//                      midDisplacement[n][m], prob_r,
//                      histDiff[n][m], prob_h,
//                      prob_t,
//                      sum_prob[n][m],
//                      coverOldByNew[n][m],coverNewByOld[n][m]
//                      );
//            }
//        }
//        
//        float velocitySum = 0;
//        int velocityBlobCount = 0;
//        float currentMaxProb = 1000;
//        int maxIterations = [newBlobs count] + [blobs count];
//        int iterations = 0;
//        NSMutableSet* usedOldBlobs = [[NSMutableSet alloc] init];
//        SPLog(@"\n-- Matching Results --\n");
//        while (iterations<=maxIterations)
//        {
//            // find the highest significant tau value
//            float maxProb = -1000;
//            int maxN=0,maxM=0;
//            for (n=0;n<[newBlobs count];n++) {
//                for (m=0;m<[blobs count];m++) {
//                    if (sum_prob[n][m]>maxProb && sum_prob[n][m]<currentMaxProb) {
//                        maxProb = sum_prob[n][m];
//                        maxN = n;
//                        maxM = m;
//                    }
//                }
//            }
//            
//            // take pick, remove paired candidates from arrays
//            Blob* oldBlob = (Blob*)[blobs objectAtIndex:maxM];
//            Blob* newBlob = (Blob*)[newBlobs objectAtIndex:maxN];
//            
//            // check mid displacement constraints and update mean velocity
//            if (midDisplacement[maxN][maxM]>maxDisplacement)
//            {
//                SPLog(@"pairing new blob #%d / old blob ID#%ld rejected. dR=%4.1f violates maximum velocity constraint.",maxN,[oldBlob tag],midDisplacement[maxN][maxM]);
//            }
//            else if (useMeanVelocityConstraint && midDisplacement[maxN][maxM] > maxMeanVelocityDisplacement)
//            {
//                SPLog(@"pairing new blob #%d / old blob ID#%ld rejected. dR=%4.1f violates mean velocity constraint.",maxN,[oldBlob tag],midDisplacement[maxN][maxM]);
//            }
//            else if ([newBlob tag]==0 && ![usedOldBlobs containsObject:oldBlob])
//            {
//                float velocity = midDisplacement[maxN][maxM]/deltaT;
//                [newBlob setVelocity:velocity];
//                if (velocity>0.5) {
//                    velocitySum+=velocity;
//                    velocityBlobCount++;
//                }
//                
//                SPLog(@"pairing new blob #%d / old blob ID=#%ld accepted, velocity %4.1f m/s",
//                      maxN,[oldBlob tag],velocity);
//                
//                [newBlob setTag:[oldBlob tag] color:[oldBlob color]];
//                [newBlob setTimestamp:[eoScan timestamp]];
//                [usedOldBlobs addObject:oldBlob];
//                [oldBlob release];
//            }
//            currentMaxProb=maxProb;
//            iterations++;
//        }
//        // update mean velocity
//        float newMeanVelocity = 0;
//        if (velocityBlobCount>0) newMeanVelocity = velocitySum/velocityBlobCount;
//        if (newMeanVelocity > 2.0) meanVelocity = newMeanVelocity;
//        
//        SPLog(@"\ncurrent mean velocity: %4.2f m/s",meanVelocity);
//        
//        // all new blobs without id now get one
//        NSMutableArray* result = [[NSMutableArray alloc] init];
//        SPLog(@"\n-- New ID assignments --");
//        for (n=0;n<[newBlobs count]; n++)
//        {
//            Blob* newBlob = [newBlobs objectAtIndex:n];
//            if ([newBlob tag]==0) {
//                [self tagBlob:newBlob];
//                SPLog(@"new blob #%d at [%3d,%3d]\t with |H|=%5d is assigned new ID#%ld",
//                      n,[newBlob midX],[newBlob midY],(int)[newBlob histogramSize],[newBlob tag]);
//            }
//            [result addObject:newBlob];
//        }
//        
//        // list those old blobs not being lined up again
//        NSEnumerator* obn = [blobs objectEnumerator];
//        Blob* ob;
//        SPLog(@"\n-- Goners --");
//        BOOL hadGoners=NO;
//        while ((ob=[obn nextObject])!=nil) {
//            if ([usedOldBlobs containsObject:ob]) continue;
//            SPLog(@"ID#%ld at [%3d,%3d]\t |H|=%5d",[ob tag],[ob midX],[ob midY],(int)[ob histogramSize]);
//            hadGoners=YES;
//        }
//        if (!hadGoners) SPLog(@"none.");
//        
//        // handle merging
//        SPLog(@"\n-- Merges --");
//        // for each new blob check coverage of old blobs
//        BOOL hadMerges = NO;
//        for (n=0;n<[result count];n++) {
//            float candidates[50];
//            int candidateCount = 0;
//            Blob* newBlob = (Blob*)[result objectAtIndex:n];
//            for (m=0;m<[blobs count];m++) {
//                if (coverOldByNew[n][m]>mergeThreshold)
//                {
//                    candidates[candidateCount++]=m;
//                }
//            }
//            // if there are more than one candidate other than self fulfilling the criteria
//            // assume a merge has happened
//            if (candidateCount>1) {
//                // remove the merged old blobs from the matching process
//                // and print them out.
//                SPLog(@"Blobs ");
//                for (m=0;m<candidateCount;m++) {
//                    SPLog(@"#%ld ",[(Blob*)[blobs objectAtIndex:candidates[m]] tag]);
//                }
//                SPLog(@" seem to have merged into blob #ID%ld, ",[newBlob tag]);
//                [self tagBlob:newBlob];
//                SPLog(@"Assigned new ID#%ld\n",[newBlob tag]);
//                hadMerges=YES;
//            }
//        }
//        if (!hadMerges) SPLog(@"none.");
//        
//        // handle splitting
//        SPLog(@"\n-- Slits --");
//        // for each new blob check coverage of old blobs
//        BOOL hadSplits = NO;
//        for (m=0;m<[blobs count];m++) {
//            float candidates[50];
//            int candidateCount = 0;
//            Blob* oldBlob = (Blob*)[blobs objectAtIndex:m];
//            for (n=0;n<[result count];n++) {
//                if (coverNewByOld[n][m]>mergeThreshold)
//                {
//                    candidates[candidateCount++]=n;
//                }
//            }
//            if (candidateCount>1) {
//                hadSplits=YES;
//                SPLog(@"Blob ID#%ld seems to have split into blobs ",[oldBlob tag]);
//                for (n=0;n<candidateCount;n++) {
//                    Blob* newBlob = (Blob*)[result objectAtIndex:candidates[n]]; 
//                    SPLog(@"#%ld ",[newBlob tag]);
//                }
//                SPLog(@"\n");
//                SPLog(@"Assigned new IDs\n");
//                for (n=0;n<candidateCount;n++) {
//                    Blob* newBlob = (Blob*)[result objectAtIndex:candidates[n]]; 
//                    unsigned long oldTag = [newBlob tag];
//                    [self tagBlob:newBlob];
//                    SPLog(@"#%ld --> #%ld\n",oldTag,[newBlob tag]);
//                }
//            }
//        }
//        if (!hadSplits) SPLog(@"none.");
//        
//        [usedOldBlobs release];
//        [newBlobs release];
//        
//        [blobs removeAllObjects];
//        [blobs addObjectsFromArray:result];
//        [blobs retain];
//        
//        // sort ascending by ID
//        // NSSortDescriptor* sorta = [[NSSortDescriptor alloc] initWithKey:@"tag" ascending:YES];
//        
//        [result sortUsingSelector: @selector(compareTag:) ];
//        //[result sortUsingDescriptors:[NSArray arrayWithObject:sorta]];
//        
//        if (doTracking) [trackStore addBlobArray:result];
//        
//        int i;
//        n=[result count];
//        for (i=0;i<n;i++) {
//            Blob* blob = (Blob*)[result objectAtIndex:i];
//            [eoScan addToBlobs:blob];
//            //	EOBlob* eoBlob = [blob eoBlob];
//            //	[eoScan addToEOBlobs:eoBlob];
//        }
//        SPLog(@"\n");
//

    }
    
};

#endif
