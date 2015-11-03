if (typeof (M3D) == 'undefined' || M3D == null) {

    M3D = {};

    /**
     * Class representing a graph of tracks constructed
     * from the data in a tracking dictionary. Also contains
     * helper methods for filtering by id and methods to
     * index or key up the data along criteria.
     */
	function TrackingGraph(dictionary, filterById, filterByTrackLength) {
	    // Data
	    this.dictionary = dictionary;
	    this.id = false;
		this.nodes = dictionary.tree.nodes;
		this.links = dictionary.tree.links;
		// Filtering and indexing
		this.linkedNodes = new Set();
		this.steps = new Array();
		this.ids = new Array();

		this.trackLengthIndex = {
		    values : new Array(),
		    map : new Map()
		}
		this.selectedTrackLength = false;

        /**
         * Sets the id to filter by.
         * @param {number} id
         */
        this.setId = function(filterById) {
            this.id = filterById;
            if (filterById) {
                this.linkedNodes = this.getLinkedNodesByID(parseInt(filterById));
            }
        }

        /**
         * Sets the track length to filter by
         * @param {number} id
         */
        this.setTrackLength = function(length) {
            var self = this;
            this.selectedTrackLength = length;
            if (length) {
                var len = parseInt(this.selectedTrackLength);
                var ids = this.trackLengthIndex.map.get(len);
                if (ids) {
                    // Add all linked nodes of all tracks of the selected length
                    ids.forEach(function(id) {
                        var li = self.getLinkedNodesByID(parseInt(id));
                        li.forEach(function(n) {
                            self.linkedNodes.add(n);
                        });
                    });
                }
            }
        }

        // Filtering

        /** @return a function to filter nodes. */
        this.nodeFilter = function(node) {
            if (this.id || this.selectedTrackLength) {
		        return this.linkedNodes.has(node);
            }
            return true;
        }

        /** @return a function to filter links. */
	    this.linkFilter = function(link) {
	        if (this.id || this.selectedTrackLength) {
                return this.linkedNodes.has(link.source) || this.linkedNodes.has(link.target);
            }
            return true;
        }

        this.filteredNodes = function() {
            var self = this;
            var fn = new Array();
            this.nodes.forEach(function(node) {
                if (self.nodeFilter(node)) {
                    fn.push(node);
                }
            });
            return fn;
        }

        this.filteredLinks = function() {
            var self = this;
            var fl = new Array();
            this.links.forEach(function(link) {
                if (self.linkFilter(link)) {
                    fl.push(link);
                }
            });
            return fl;
        }

        // Indexing and Searching

        /**
         * Turns the uuid references in the links array into
         * object references.
         * @param {array} nodes
         * @param {array} links
         */
        this.resolveLinks = function() {
            var self = this;
            var links = this.links;
            this.nodes.forEach(function(node) {
                self.links.forEach(function(link) {
                    if (link.source == node.uuid) {
                        link.source = node;
                    }
                    if (link.target == node.uuid) {
                        link.target = node;
                    }
                });
            });
        };

        /**
         * Gets all nodes with the matching id
         * @param {Array} nodes
         * @param {number} id
         * @return {Set} of nodes.
         */
        this.nodesWithId = function(id) {
            var nodeSet = new Set();
            this.nodes.forEach(function(node) {
                if (node.id == id)
                    nodeSet.add(node);
            });
            return nodeSet;
        }

        /**
         * Obtain by id a list of nodes that are linked, starting
         * with all nodes with the given id.
         * @param {Array} nodes
         * @param {Array} links
         * @param {number} id
         * @return {Set} of nodes.
         */
        this.getLinkedNodesByID = function(id) {
            var self = this;
            var roots = this.nodesWithId(id);
            var linked = new Set();
            roots.forEach(function(root) {
                self.getLinkedNodes(root, linked);
            });
            return linked;
        }

        /**
         * Recursively collects all ids of nodes linking
         * to or linked to any node with this id.
         */
        this.getLinkedNodes = function(node,linkedNodes) {
            var self = this;
            if (linkedNodes.has(node)) {
                return;
            }
            linkedNodes.add(node);
            var nextNodes = this.nodesWithId(node.id);
            // Now gather up all nodes that point to or
            // from any of the nodes in the list
            this.links.forEach(function(link) {
                if (nextNodes.has(link.target) || nextNodes.has(link.source)) {
                    var nextNode = null;
                    if (link.target == node) {
                        self.getLinkedNodes(link.source, linkedNodes);
                    } else if (link.source == node) {
                        self.getLinkedNodes(link.target, linkedNodes);
                    }
                }
            });
        }

        /**
         * Indexes the whole tracking dictionary into a few
         * indexes used for navigation.
         */
        this.index = function() {
            var self = this;
            // Index ids/steps
            var ids = this.ids;
            self.nodes.forEach(function(node) {
                if (self.steps.indexOf(node.step) < 0)
                    self.steps.push(node.step);
                if (self.ids.indexOf(node.id) < 0)
                    self.ids.push(node.id);
                });
            this.ids.sort(function(a,b) {
                return parseInt(a.id) < parseInt(b.id);
            });
            this.steps.sort(function(a,b) {
                return parseInt(a.step) < parseInt(b.step);
            });
            // index track dictionary
            this.dictionary.tracks.forEach(function(track) {
                if (self.trackLengthIndex.values.indexOf(track.length) < 0) {
                    self.trackLengthIndex.values.push(track.length);
                }
                var trackIds = self.trackLengthIndex.map.get(track.length);
                if (!trackIds) {
                    trackIds = new Set();
                    self.trackLengthIndex.map.set(track.length,trackIds);
                }
                trackIds.add(track.id);
            });
            // Sort descending
            this.trackLengthIndex.values.sort(function(a,b) {
                return (a > b);
            });
            console.log(this.trackLengthIndex.values);
        };

        // Finally index the data
        this.resolveLinks();
		this.index();
		this.setId(filterById);
		this.setTrackLength(filterByTrackLength)
    }
    M3D.TrackingGraph = TrackingGraph;

    /**
     * Courtesy of https://css-tricks.com/snippets/javascript/get-url-variables/
     */
    M3D.getUrlParameter= function(param)
    {
       var query = window.location.search.substring(1);
       var vars = query.split("&");
       for (var i=0;i<vars.length;i++) {
               var pair = vars[i].split("=");
               if(pair[0] == param){return pair[1];}
       }
       return(false);
    }

    /**
     * Reloads the current page with an id for filtering.
     * @param A node or link or any other object with property 'id'.
     */
    M3D.filterById = function(d) {
        var location = window.location.toString();
        var i = location.indexOf("?");
        if ( i > 0) {
            location = location.substring(0,i);
        }
        window.location = location + "?id=" + d;
    }

    M3D.filterByTrackLength = function(d) {
        var location = window.location.toString();
        var i = location.indexOf("?");
        if ( i > 0) {
            location = location.substring(0,i);
        }
        window.location = location + "?tracklength=" + d;
    }
}