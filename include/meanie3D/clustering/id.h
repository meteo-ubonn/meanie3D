#ifndef M3D_ID
#define M3D_ID

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <map>
#include <set>
#include <vector>

extern "C"
{
    namespace m3D {

        /** Data type for object identifiers
         */
        typedef unsigned long           id_t;

        typedef std::set<id_t>          id_set_t;

        typedef std::vector<id_t>       id_vec_t;

        typedef std::map<id_t,id_set_t> id_map_t;

        // Constants

        static const id_t NO_ID = std::numeric_limits<id_t>::max();

        static const id_t MIN_ID = 0;

        static const id_t MAX_ID = std::numeric_limits<id_t>::max() - 1;

        // Methods

        /** Increment the ID, rotating around if necessary.
         * @param current id, incremented after the call
         * @return next id
         */
        m3D::id_t next_id(m3D::id_t &current)
        {
            if (current >= m3D::MAX_ID)
            {
                current = m3D::MIN_ID;
            }
            else
            {
                current++;
            }

            return current;
        }
    }
}

#endif
