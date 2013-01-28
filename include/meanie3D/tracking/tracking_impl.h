#ifndef _M3D_TRACKING_CLASS_Impl_H_
#define _M3D_TRACKING_CLASS_Impl_H_

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <vector>
#include <netcdf>
#include <set>

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
    Tracking<T>::track(typename ClusterList<T>::ptr previous,
                       typename ClusterList<T>::ptr current,
                       const NcVar &track_variable,
                       Verbosity verbosity)
    {
        // sanity check
        
        if ( current->clusters.size() == 0 )
        {
            if ( verbosity >= VerbosityNormal )
            {
                cout << "No current clusters. Skipping tracking." << endl;
                
                return;
            }
        }
        
        
        // figure out tracking variable index
        
        int index = index_of_first( current->feature_variables, track_variable );
        
        assert( index >= 0 );
        
        size_t tracking_var_index = (size_t) index;
        
        
        // Valid min/max
        
        T valid_min, valid_max;
        
        track_variable.getAtt("valid_min").getValues( &valid_min );

        track_variable.getAtt("valid_max").getValues( &valid_max );
        
        
        // Find the highest id 
        
        size_t current_id = 0;

        // Find out the highest id from the previous list
        
        typename Cluster<T>::list::iterator ci;
        
        for ( ci = previous->clusters.begin(); ci != previous->clusters.end(); ci++ )
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
        // TODO: add timestamp to the cluster file format
        // TODO: timestamp in the original file other than in the filename? Spec?
        
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

        T maxDisplacement = this->m_maxVelocity * this->m_deltaT;
        
        if ( verbosity >= VerbosityDetails )
        {
            printf("max velocity  constraint at %4.1f m/s at deltaT %4.0fs -> dR_max = %7.1fm\n",
                   this->m_maxVelocity, this->m_deltaT, maxDisplacement );
        }
        
//        T maxMeanVelocityDisplacement = 0.0;
//        
//        if (m_useMeanVelocityConstraint)
//        {
//            m_useMeanVelocityConstraint = deltaT * meanVelocity * m_meanVelocitySecurityPercentage;
//            
//            printf("mean velocity constraint at %4.2f m/s, dR_max = %7.1f",
//                   meanVelocity * m_meanVelocitySecurityPercentage, maxMeanVelocityDisplacement );
//        }
        
//        SPLog(@"\n");
        
        // Prepare the newcomers for re-identification
        current->erase_identifiers();

        // Get the counts
        
        size_t old_count = previous->clusters.size();
        size_t new_count = current->clusters.size();
        
        // Create result matrices
        
        matrix_t rank_correlation = create_matrix(new_count,old_count);
        matrix_t midDisplacement = create_matrix(new_count,old_count);
        matrix_t histDiff = create_matrix(new_count,old_count);
        matrix_t sum_prob = create_matrix(new_count,old_count);
        matrix_t coverOldByNew = create_matrix(new_count,old_count);
        matrix_t coverNewByOld = create_matrix(new_count,old_count);
        
        size_t maxHistD = numeric_limits<size_t>::min();
        
        T maxMidD = numeric_limits<T>::min();
        
        // compute correlation matrix
        
        int n,m;
        
        // check out the values for midDisplacement, histSizeDifference and kendall's tau
        
        for ( n=0; n < new_count; n++ )
        {
            typename Cluster<T>::ptr newCluster = current->clusters[n];
            
            typename Histogram<T>::ptr newHistogram = newCluster->histogram(tracking_var_index,valid_min,valid_max);
            
            for ( m=0; m < old_count; m++ )
            {
                typename Cluster<T>::ptr oldCluster = previous->clusters[m];
                
                typename Histogram<T>::ptr oldHistogram = oldCluster->histogram(tracking_var_index,valid_min,valid_max);
                
                // calculate spearman's tau
                
                rank_correlation[n][m] = newHistogram->correlate_kendall( oldHistogram );
                
                // calculate average mid displacement
                
                vector<T> dx = newCluster->weighed_center(current->spatial_dimension,tracking_var_index)
                             - oldCluster->weighed_center(current->spatial_dimension,tracking_var_index);

                midDisplacement[n][m] = vector_norm(dx);

                size_t max_size = max( newHistogram->sum(), oldHistogram->sum() );
                
                histDiff[n][m] = (max_size==0) ? 1.0 : (abs( (T)newHistogram->sum() - (T)oldHistogram->sum() ) / ((T)max_size) );
                
                coverNewByOld[n][m] = newCluster->percent_covered_by( oldCluster );
                
                coverOldByNew[n][m] = oldCluster->percent_covered_by( newCluster );
                
                // track maxHistD and maxMidD
                
                if ( histDiff[n][m] > maxHistD )
                {
                    maxHistD = histDiff[n][m];
                }
                
                if (midDisplacement[n][m]>maxMidD)
                {
                    maxMidD = midDisplacement[n][m];
                }

            }
        }
        
        // Can't have zeros here
        
        if (maxHistD==0) maxHistD = 1;
        
        // normalise maxMidD with maximum possible distance
        // why?
        
        // maxMidD = sqrt( 2*200*200 );

        // calculate the final probability
        cout << endl << "-- Correlation Table --" << endl;
        
        for ( n=0; n < new_count; n++ )
        {
            typename Cluster<T>::ptr newCluster = current->clusters[n];
            
            if ( verbosity >= VerbosityDetails )
            {
                cout << endl << "Correlating new Blob #" << n
                << " (histogram size " << newCluster->histogram(tracking_var_index,valid_min,valid_max)->sum() << ")"
                << " at " << newCluster->mode << endl;
            }

            for ( m=0; m < old_count; m++ )
            {
                typename Cluster<T>::ptr oldCluster = previous->clusters[m];

                float prob_r = m_dist_weight * erfc( midDisplacement[n][m] / maxMidD );
                float prob_h = m_size_weight * erfc( histDiff[n][m] / maxHistD );
                float prob_t = m_corr_weight * rank_correlation[n][m];
                sum_prob[n][m] = prob_t + prob_r + prob_h;
                
                if ( verbosity >= VerbosityDetails )
                {
                    printf("\t<ID#%4llu>:\t(|H|=%5lu)\t\tdR=%4.1f (%5.4f)\t\tdH=%5.4f (%5.4f)\t\ttau=%7.4f\t\tsum=%6.4f\t\t\tcovON=%3.2f\t\tcovNO=%3.2f\n",
                           oldCluster->id,
                           oldCluster->histogram(tracking_var_index,valid_min,valid_max)->sum(),
                           midDisplacement[n][m],
                           prob_r,
                           histDiff[n][m],
                           prob_h,
                           rank_correlation[n][m],
                           sum_prob[n][m],
                           coverOldByNew[n][m],
                           coverNewByOld[n][m]);
                }
            }
        }
        
        float velocitySum = 0;
        
        int velocityBlobCount = 0;
        
        float currentMaxProb = numeric_limits<float>::max();
        
        int maxIterations = new_count * old_count;
        
        int iterations = 0;
        
        // NSMutableSet* usedOldBlobs = [[NSMutableSet alloc] init];
        
        set< typename Cluster<T>::ptr > used_clusters;
        
        printf("\n-- Matching Results --\n");
        
        while ( iterations <= maxIterations )
        {
            // find the highest significant correlation match
            
            float maxProb = numeric_limits<float>::min();
            
            size_t maxN=0, maxM=0;
            
            for ( n=0; n < new_count; n++ )
            {
                for ( m=0; m < old_count; m++ )
                {
                    if ( (sum_prob[n][m] > maxProb) && (sum_prob[n][m] < currentMaxProb) )
                    {
                        maxProb = sum_prob[n][m];

                        maxN = n;
                        
                        maxM = m;
                    }
                }
            }

            // take pick, remove paired candidates from arrays
            
            typename Cluster<T>::ptr new_cluster = current->clusters[maxN];

            typename Cluster<T>::ptr old_cluster = previous->clusters[maxM];
            
            // check mid displacement constraints and update mean velocity
            
            if ( midDisplacement[maxN][maxM] > maxDisplacement )
            {
                if ( verbosity >= VerbosityDetails )
                {
                    printf("pairing new blob #%4d / old blob ID#%lld rejected. dR=%4.1f violates maximum velocity constraint.\n", maxN, old_cluster->id, midDisplacement[maxN][maxM] );
                }
            }
//            else if (m_useMeanVelocityConstraint && midDisplacement[maxN][maxM] > maxMeanVelocityDisplacement )
//            {
//                if ( verbosity >= VerbosityDetails )
//                {
//                    printf("pairing new blob #%d / old blob ID#%ld rejected. dR=%4.1f violates mean velocity constraint.",
//                           maxN, matchedPrev->id, midDisplacement[maxN][maxM] );
//                }
//            }
            else if ( new_cluster->id == Cluster<T>::NO_ID )
            {
                // new cluster not tagged yet
                
                if ( used_clusters.find(old_cluster) == used_clusters.end() )
                {
                    // old cluster not matched yet

                    float velocity = midDisplacement[maxN][maxM] / this->m_deltaT;
                
                    velocitySum += velocity;
                        
                    velocityBlobCount++;
                    
    //                if ( verbosity >= VerbosityDetails )
    //                {
                        printf("pairing new blob #%4lu / old blob ID=#%4llu accepted, velocity %4.1f m/s\n", maxN, old_cluster->id, velocity );
    //                }
                    
                    new_cluster->id = old_cluster->id;
                    
                    used_clusters.insert( old_cluster );
                    
                    current->tracked_ids.push_back( new_cluster->id );
                }
            }
            
            currentMaxProb = maxProb;
            
            iterations++;
        }

//        // update mean velocity
//        float newMeanVelocity = 0;
//        if (velocityBlobCount>0) newMeanVelocity = velocitySum/velocityBlobCount;
//        if (newMeanVelocity > 2.0) meanVelocity = newMeanVelocity;
//
//        SPLog(@"\ncurrent mean velocity: %4.2f m/s",meanVelocity);
        
        // all new blobs without id now get one
        
        printf("\n-- New ID assignments --\n");
        
        for (n=0;n<new_count;n++)
        {
            typename Cluster<T>::ptr c = current->clusters[n];
            
            if ( c->id == Cluster<T>::NO_ID )
            {
                c->id = current_id;
                
                current->new_ids.push_back( c->id );
                
                cout << "new cluster #" << n << " at " << c->mode << " with |H|=" << c->histogram(tracking_var_index,valid_min,valid_max)->sum()
                     << " is assigned new ID#" << c->id << endl;
                
                current_id++;
            }
        }

        // list those old blobs not being lined up again
        bool hadGoners = false;
        cout << "\n-- Goners --" << endl;
        for ( ci = previous->clusters.begin(); ci != previous->clusters.end(); ci++ )
        {
            typename Cluster<T>::ptr c = *ci;
            
            typename set< typename Cluster<T>::ptr >::iterator f = used_clusters.find( c );
            
            if ( f != used_clusters.end() ) continue;
            
            cout << "ID#" << c->id << " at " << c->mode << " |H|=" << c->histogram(tracking_var_index,valid_min,valid_max)->sum() << endl;
            
            current->dropped_ids.push_back(c->id);

            hadGoners = true;
            
        }
        if (!hadGoners)
        {
            cout << "none." << endl;
        }
        
//        NSEnumerator* obn = [blobs objectEnumerator];
//        Blob* ob;
//        SPLog(@"\n-- Goners --");
//        BOOL hadGoners=NO;
//        while ((ob=[obn nextObject])!=nil)
//        {
//            if ([usedOldBlobs containsObject:ob]) continue;
//            SPLog(@"ID#%ld at [%3d,%3d]\t |H|=%5d",[ob tag],[ob midX],[ob midY],(int)[ob histogramSize]);
//            hadGoners=YES;
//        }
//        if (!hadGoners) SPLog(@"none.");
        
        
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
        
        current->tracking_performed = true;
    }
    
};

#endif
