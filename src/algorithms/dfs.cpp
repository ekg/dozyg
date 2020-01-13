#include "dfs.hpp"

namespace gyeet {
namespace algorithms {

using namespace handlegraph;


// depth first search across node traversals with interface to traversal tree via callback
void dfs(
    const HandleGraph& graph,
    const std::function<void(const handle_t&)>& handle_begin_fn,  // called when node orientation is first encountered
    const std::function<void(const handle_t&)>& handle_end_fn,    // called when node orientation goes out of scope
    const std::function<bool(const handle_t&)>& handle_skip_fn,   // tells us if we should skip the handle
    const std::function<bool(void)>& break_fn,                    // called to check if we should stop the DFS; we stop when true is returned.
    const std::function<void(const edge_t&)>& edge_fn,            // called when an edge is encountered
    const std::function<void(const edge_t&)>& tree_fn,            // called when an edge forms part of the DFS spanning tree
    const std::function<void(const edge_t&)>& edge_curr_fn,       // called when we meet an edge in the current tree component
    const std::function<void(const edge_t&)>& edge_cross_fn,      // called when we meet an edge in an already-traversed tree component
    const std::vector<handle_t>& sources,                         // start only at these node traversals
    const ska::flat_hash_set<handle_t>& sinks                     // when hitting a sink, don't keep walking
    ) {

    // to maintain search state
    enum SearchState { PRE = 0, CURR, POST };
    ska::flat_hash_map<handle_t, SearchState> state; // implicitly constructed entries will be PRE.

    // to maintain stack frames
    struct Frame {
        handle_t handle;
        std::vector<edge_t>::iterator begin, end;
        Frame(handle_t h,
              std::vector<edge_t>::iterator b,
              std::vector<edge_t>::iterator e)
            : handle(h), begin(b), end(e) { }
    };

    // maintains edges while the node traversal's frame is on the stack
    ska::flat_hash_map<handle_t, std::vector<edge_t> > edges;

    // do dfs from given root.  returns true if terminated via break condition, false otherwise
    std::function<bool(const handle_t&)> dfs_single_source = [&](const handle_t& root) {
        
#ifdef debug
        cerr << "Starting a DFS from " << graph.get_id(root) << " " << graph.get_is_reverse(root) << endl;
#endif
        
        // to store the stack frames
        std::deque<Frame> todo;
        if (state[root] == SearchState::PRE) {
#ifdef debug
            cerr << "\tMoving from PRE to CURR" << endl;
#endif
        
            state[root] = SearchState::CURR;
                
            // Collect all the edges attached to the outgoing side of the
            // traversal.
            auto& es = edges[root];
            // follow edges?
            graph.follow_edges(
                root, false,
                [&](const handle_t& next) {
                    if (!handle_skip_fn(next)) {
                        es.push_back(graph.edge_handle(root, next));
#ifdef debug
                        cerr << "\t\tWill need to investigate edge to "
                             << graph.get_id(next) << " " << graph.get_is_reverse(next) << endl;
#endif
                    }
                });
                
            todo.push_back(Frame(root, es.begin(), es.end()));
            // run our discovery-time callback
            handle_begin_fn(root);
            // and check if we should break
            if (break_fn()) {
                return false;
            }
        }
        // now begin the search rooted at this NodeTraversal
        while (!todo.empty()) {
            // get the frame
            auto& frame = todo.back();
            // and set up reference to it
            auto handle = frame.handle;
            auto edges_begin = frame.begin;
            auto edges_end = frame.end;
            todo.pop_back();
            
#ifdef debug
            cerr << "Continuing DFS; at " << graph.get_id(handle) << " " << graph.get_is_reverse(handle) << endl;
#endif
            
            // run through the edges to handle
            while (edges_begin != edges_end) {
                auto& edge = *edges_begin;
                // run the edge callback
                edge_fn(edge);
                    
                // what's the handle we'd get to following this edge
                const handle_t& target = graph.traverse_edge_handle(edge, handle);
                
#ifdef debug
                cerr << "\tProcessing edge " << graph.get_id(edge.first) << " " << graph.get_is_reverse(edge.first) << " to "
                    << graph.get_id(edge.second) << " " << graph.get_is_reverse(edge.second) << endl;
                cerr << "\tHandled as edge " << graph.get_id(handle) << " " << graph.get_is_reverse(handle) << " to "
                    << graph.get_id(target) << " " << graph.get_is_reverse(target) << endl;
#endif

                auto search_state = state[target];
                // if we've not seen it, follow it
                if (search_state == SearchState::PRE) {

#ifdef debug
                    cerr << "\t\tDiscovered new destination!" << endl;
#endif
                
                    tree_fn(edge);
                    // save the rest of the search for this handle on the stack
                    todo.push_back(Frame(handle, ++edges_begin, edges_end));
                    
#ifdef debug
                    cerr << "\t\tPausing existing stack frame (" << (edges_end - edges_begin) << " edges left)" << endl;
#endif
                    
                    // switch our focus to the handle at the other end of the edge
                    handle = target;
                    
#ifdef debug
                    cerr << "\t\tFocusing and moving " << graph.get_id(handle) << " "
                        << graph.get_is_reverse(handle) << " from PRE to CURR" << endl;
#endif
                    
                    // and store it on the stack
                    state[handle] = SearchState::CURR;
                    auto& es = edges[handle];

                    // only walk out of handle that are not the sink
                    if (sinks.empty() || sinks.count(handle) == false) {
                        // follow edges?
                        graph.follow_edges(
                            handle, false,
                            [&](const handle_t& next) {
                                if (!handle_skip_fn(next)) {
                                    es.push_back(graph.edge_handle(handle, next));
#ifdef debug
                                    cerr << "\t\t\tWill need to investigate edge to "
                                         << graph.get_id(next) << " " << graph.get_is_reverse(next) << endl;
#endif
                                }
                            });
                    }
                    edges_begin = es.begin();
                    edges_end = es.end();
                    // run our discovery-time callback
                    handle_begin_fn(handle);

                    if (break_fn()) {
                        return false;
                    }
                    
                    // We will continue with edges_begin and edges_end and handle adjusted to this new loaded stack frame.
                    
#ifdef debug
                    cerr << "Now continuing DFS from " << graph.get_id(handle) << " " << graph.get_is_reverse(handle) << endl;
#endif
                    
                } else if (search_state == SearchState::CURR) {
                    // if it's on the stack
                    
#ifdef debug
                    cerr << "\tDestination already discovered" << endl;
#endif
                    
                    edge_curr_fn(edge);
                    ++edges_begin;
                } else {
                    // it's already been handled, so in another part of the tree
                    
#ifdef debug
                    cerr << "\tDestination already completed by a shorter path" << endl;
#endif
                    
                    edge_cross_fn(edge);
                    ++edges_begin;
                }
            }
            state[handle] = SearchState::POST;
            handle_end_fn(handle);
            edges.erase(handle); // clean up edge cache
        }

        return false;
    };

    if (sources.empty()) {
        // attempt the search rooted at all NodeTraversals
        graph.for_each_handle([&](const handle_t& handle_fwd) {
                handle_t handle_rev = graph.flip(handle_fwd);
                dfs_single_source(handle_fwd);
                dfs_single_source(handle_rev);
            });
    } else {
        for (auto source : sources) {
            dfs_single_source(source);
        }
    }
}

void dfs(const HandleGraph& graph,
         const std::function<void(const handle_t&)>& handle_begin_fn,
         const std::function<void(const handle_t&)>& handle_end_fn,
         const std::vector<handle_t>& sources,
         const ska::flat_hash_set<handle_t>& sinks) {
    auto edge_noop = [](const edge_t& e) { };
    dfs(graph,
        handle_begin_fn,
        handle_end_fn,
        [](const handle_t& h) { return false; },
        [](void) { return false; },
        edge_noop,
        edge_noop,
        edge_noop,
        edge_noop,
        sources,
        sinks);
}

void dfs(const HandleGraph& graph,
         const std::function<void(const handle_t&)>& handle_begin_fn,
         const std::function<void(const handle_t&)>& handle_end_fn,
         const std::function<bool(const handle_t&)>& handle_skip_fn,
         const std::function<bool(void)>& break_fn,
         const std::vector<handle_t>& sources) {
    auto edge_noop = [](const edge_t& e) { };
    ska::flat_hash_set<handle_t> empty_sinks;
    dfs(graph,
        handle_begin_fn,
        handle_end_fn,
        handle_skip_fn,
        break_fn,
        edge_noop,
        edge_noop,
        edge_noop,
        edge_noop,
        sources,
        empty_sinks);
}

void dfs(const HandleGraph& graph,
         const std::function<void(const handle_t&)>& handle_begin_fn,
         const std::function<void(const handle_t&)>& handle_end_fn,
         const std::function<bool(void)>& break_fn) {
    auto edge_noop = [](const edge_t& e) { };
    std::vector<handle_t> empty_sources;
    ska::flat_hash_set<handle_t> empty_sinks;
    dfs(graph,
        handle_begin_fn,
        handle_end_fn,
        [](const handle_t& h) { return false; },
        break_fn,
        edge_noop,
        edge_noop,
        edge_noop,
        edge_noop,
        empty_sources,
        empty_sinks);
}

}
}
