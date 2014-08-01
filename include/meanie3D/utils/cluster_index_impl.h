#ifndef _M3D_ClusterIndex_Impl_H_
#define _M3D_ClusterIndex_Impl_H_

namespace m3D { namespace utils {

    template <typename T>
    ClusterIndex<T>::ClusterIndex(typename Cluster<T>::list &list,
                                  const vector<size_t> &dimensions)
    {
        // Figure out the rank of the multiarray needed
        // from the first cluster's mode
        
        this->m_index = new MultiArrayBlitz<cfa::id_t>(dimensions);
        
        if (!list.empty())
        {
            for (size_t ci=0; ci < list.size(); ci++)
            {
                typename Cluster<T>::ptr c = list.at(ci);
                
                for (size_t pi=0; pi < c->points.size(); pi++)
                {
                    typename Point<T>::ptr p = c->points.at(pi);
                    
                    this->m_index->set(p->gridpoint,c->id);
                }
            }
        }
    }
    
    template <typename T>
    ClusterIndex<T>::~ClusterIndex()
    {
        if (m_index != NULL)
        {
            delete m_index;
            m_index = NULL;
        }
    }
    
    template <typename T>
    typename ClusterIndex<T>::index_t *
    ClusterIndex<T>::data()
    {
        return this->m_index;
    }
    
    
    template <typename T>
    size_t
    ClusterIndex<T>::count_common_points(const ClusterIndex<T> *a,
                                         const ClusterIndex<T> *b,
                                         ::cfa::id_t id)
    {
        assert(a->data()->rank() == b->data()->rank());
        
        size_t common_points;
        
        size_t dim = a->data().rank();

        vector<int> gridpoint(dim);
        
        if (dim==2)
        {
            for (size_t i1=0; i1 < a->data()->get_dimensions(0); i1++)
            {
                gridpoint[0] = i1;

                for (size_t i2=0; i2 < a->data()->get_dimensions(1); i2++)
                {
                    gridpoint[1] = i2;
                    
                    cfa::id_t id_a = a->data()->get(gridpoint);
                    
                    if (id_a == id)
                    {
                        cfa::id_t id_b = b->data()->get(gridpoint);
                        
                        if (id_b == id_a) common_points++;
                    }
                }
            }
        }
        else
        {
            throw std::runtime_error("Not implemented");
        }
        
        return common_points;
    }
    
    template <typename T>
    double
    ClusterIndex<T>::occupation_ratio(const Cluster<T> *cluster_a,
                                      const Cluster<T> *cluster_b) const
    {
        int common_points = 0;
        
        for (int i=0; i<cluster_a->points.size(); i++)
        {
            typename Point<T>::ptr p = cluster_a->points.at(i);
            
            if (m_index->get(p->gridpoint) == cluster_b->id)
            {
                common_points++;
            }
        }
        
        double ratio = ((double)common_points)/((double)cluster_a->points.size());
        
        return ratio;
    }

    
}}

#endif
