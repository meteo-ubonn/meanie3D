if (typeof (M3D) == 'undefined' || M3D == null) {

    M3D = {};

    /**
     * Class representing a graph of tracks constructed
     * from the data in a tracking dictionary. Also contains
     * helper methods for filtering by id and methods to
     * index or key up the data along criteria.
     */
	function TrackingGraph(dictionary, filterById) {
	    // Data
	    this.dictionary = dictionary;
	    this.id = false;
		this.nodes = dictionary.tree.nodes;
		this.links = dictionary.tree.links;
		// Filtering and indexing
		this.linkedNodes = new Set();
		this.steps = new Array();
		this.ids = new Array();
		this.tracksByLength = new Map();
		this.clustersBySize = new Map();
		this.clustersByFilenane = new Map();

        /**
         * Sets the id to filter by.
         * @param {number} id
         */
        this.setId = function(filterById) {
            this.id = filterById;
            this.linkedNodes = new Set();
            if (filterById) {
                this.linkedNodes = this.getLinkedNodesByID(parseInt(filterById));
            }
        }

        // Filtering

        /** @return a function to filter nodes. */
        this.nodeFilter = function(node) {
		    return this.id ? this.linkedNodes.has(node) : true;
        }

        /** @return a function to filter links. */
	    this.linkFilter = function(link) {
            return this.id ? this.linkedNodes.has(link.source)
            || this.linkedNodes.has(link.target) : true;
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
         * Indexes the whole tracking dictionary into a few
         * indexes used for navigation.
         */
        this.index = function() {
            var self = this;
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

        // Finally index the data
        this.resolveLinks();
		this.index();
		this.setId(filterById);
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
    M3D.navigateTo = function(d) {
        var location = window.location.toString();
        var i = location.indexOf("?");
        if ( i > 0) {
            location = location.substring(0,i);
        }
        window.location = location + "?id=" + d.id;
    }

    M3D.resizeSvg = function () {
        var width = Math.max(document.documentElement.clientWidth, window.innerWidth || 0)
        var height = Math.max(document.documentElement.clientHeight, window.innerHeight || 0)
        d3.select("svg")
            .attr("width", width)
            .attr("height", height);
    }
}