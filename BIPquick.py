import numpy as np
np.set_printoptions(threshold=np.inf, linewidth=np.inf)
import csv
import sys
import FileSettings


def main():
    if __name__ == "__main__":
        print("Entering Main()")
    sys.setrecursionlimit(60000)
    """
    The main function puts the tools from the other sections together into a cogent workflow
    """
    global nodes, links, visit_network_mask, visited_flag, subnetwork_container_nodes, subnetwork_container_links, \
        processors, counter, partition_boolean, partition_threshold, nrows, lrows
    # The .csv files containing link-node information abstracted from the SWMM file are opened and parsed into the nodes
    # and links array.  The starting point for the phantom node/link's names (should they be necessary) is determined.
    nodes_info = FileSettings.settingsdict['nodes_info']
    links_info = FileSettings.settingsdict['links_info']
    array_maker(nodes_info, links_info)
    #print(nodes)
    #print(links)
    phantom_naming_convention()

    # The visit_network_mask, unlike the other boolean list, is only initialized once.  This ensures that
    # each link can only exist on 1 subnetwork.
    visit_network_mask = [False for node in nodes]
    partition_boolean = [False for link in links]

    # The nr_directweight_u column is updated from default values to values determined by the weight contributed to each
    # node from directly upstream.  The nr_totalweight_u column is also initialized to determine max_weight, used in the
    # calculation of the static partition_threshold.
    local_node_weighting()
    null_value_convert(nodes[:, nr_totalweight_u])

    #nr_totalweight_assigner()
    print(nodes)
    print(links)


    # Using the number of processors available, the containers for the subnetwork information are initialized
    processors = FileSettings.settingsdict['multiprocessors']
    # And the static variable partition_threshold is calculated
    processors_counter = 0
    # The counter is the suffix for the phantom elements to ensure they are all uniquely named
    counter = 1
    subnetwork_container_nodes = np.empty((processors, len(nodes), len(nodes[1])), dtype=object)
    subnetwork_container_links = np.empty((processors, len(links), len(links[1])), dtype=object)

    # This loop populates the subnetwork_container_nodes depth-wise.  This is primarily for spot-checking purposes; when
    # the containers are initialized as processors x rows x columns then it's easy to track when printed.
    for i in range(processors):
        # The cumulative weight for each node needs to be reset after each subnetwork is taken out
        nodes[:, nr_totalweight_u] = 0.0
        # Determine the new cumulative weight for each node for a given network
        nr_totalweight_assigner()
        print(nodes[:, nr_directweight_u])
        print(nodes[:, nr_totalweight_u])
        partition_threshold = max_weight / (processors - processors_counter)

        # The effective_root is the root where the partition_threshold is found exactly.  Alternatively, if the
        # partition_threshold isn't found at a node or spanned by a link, then the effective_root is the nearest
        # overestimate of the partition_threshold by the cumulative node weight.
        effective_root = ideal_partition_check(processors_counter)
        # If the partition_threshold exists at a node, separate out the upstream network (populating the
        # subnetwork_container_nodes), and rejoice.
        if ideal_exists:
            print("Ideal! " + effective_root)
            #print("Entering subnetwork_carving()")
            subnetwork_carving(effective_root, i)
        # Otherwise, check to see if the partition_threshold is spanned by any link
        else:
            spanning_link = spanning_check(partition_threshold)
            print(spanning_link)
            # If the spanning_link is NOT empty (i.e. the partition_threshold is spanned) determine the location along
            # the link that the partition_threshold must exit, and put a phantom node there.  This also has the effect
            # of shortening the spanning link to upstream of the phantom node, and creating a phantom link to connect
            # the phantom node to the rest of the network.
            while spanning_link == '':
                print('Not Ideal and Not Spanned after ' + str(i) + ' partitions.')
                """
                effective_root - the nearest overestimate of the partition threshold
                Process:
                - add the effective_root to the nodes container but don't turn off its visited flag
                - find an upstream neighbor of the effective_root, save the nr_totalweight_u for the upstream neighbor
                - call subnetwork_carving(upstream_neighbor)
                - reduce the partition threshold by the nr_totalweight_u from the upstream_neighbor
                - search for spanning_links again... theres going to have to be a while loop in here somewhere
                """
                print("The effective root is: " + effective_root)
                for nrow in range(len(nodes)):
                    if nodes[nrow, ni_idx] == effective_root:
                        subnetwork_container_nodes[i, nrow] = nodes[nrow]
                        for lrow in range(len(links)):
                            if links[lrow, li_Mnode_d] == effective_root:
                                upstream_link = links[lrow, li_idx]
                                upstream_link_length = links[lrow, lr_Length]
                                upstream_node = links[lrow, li_Mnode_u]
                                for nrows2 in range(len(nodes)):
                                    if nodes[nrows2, ni_idx] == upstream_node:
                                        upstream_weight = nodes[nrows2, nr_totalweight_u]
                                        total_clipped_weight = float(upstream_weight) + weighting_function(
                                            float(upstream_link_length))
                                        break
                                break
                        nodes[nrow, nr_directweight_u] = nodes[nrow, nr_directweight_u] - total_clipped_weight
                        break
                partition_threshold = partition_threshold - total_clipped_weight
                print("The amount of weight removed from the system is: " + str(total_clipped_weight))
                print("The new partition_threshold is: " + str(partition_threshold))
                subnetwork_carving(upstream_node, i)
                spanning_link = spanning_check(partition_threshold)
                print("The spanning link is:" + spanning_link)

            length_from_start = linear_interpolator(partition_threshold, spanning_link)
            phantom_node_generator(spanning_link, length_from_start, i)
            # Also, separate out the upstream subnetwork from the phantom node, and add it to the nodes container.
            #print("Entering subnetwork_carving()")
            subnetwork_carving(nodes[phantom_array_loc, ni_idx], i)
            #print(subnetwork_container_nodes)
            print(visit_network_mask)

        for row in range(len(visit_network_mask)):
            # if nodes[row, ni_idx] in subnetwork_container_nodes[i, :, ni_idx]:
            #     nodes[row, nr_directweight_u] = 0
            if visit_network_mask[row] is True:
                nodes[row, nr_directweight_u] = 0


        print(nodes)
        print(links)

        # Update the processors_counter, which is used to determine the partition_threshold
        processors_counter = processors_counter + 1
        #quit()

    return


""" 
Globals Initialization:

This section initializes the headers that are used in the two imported arrays.  These globals are used to
specify the columns of the nodes and links array masks that are used in subsequent functions
"""
global nodes, node_id, ni_idx, ni_node_type, ni_node_type, ni_N_link_u, ni_N_link_d, ni_Mlink_u1, ni_Mlink_u2, \
    ni_Mlink_u3, ni_Mlink_d1, ni_Mlink_d2, ni_Mlink_d3, nr_directweight_u, nr_totalweight_u, link_id, li_idx, \
    li_link_type, li_geometry, li_Mnode_u, li_Mnode_d, lr_Length, lr_Slope, lr_Roughness, lr_InitialFlowrate, \
    lr_InitialUpstreamDepth, lr_InitialDnstreamDepth, lr_Target

global subnetwork_container_nodes, visit_network_mask, visited_flag, max_weight, partition_threshold, ideal_exists, \
   subnetwork_container_links, visited_flag_weight, weight_index, nodes, links


nleft = 0
node_id = 1  # node ID
ni_idx = 2  # the node index (i.e. its position in the array)
ni_node_type = 3  # int representation of the type of node
ni_N_link_u = 4  # the number of links directly upstream of this node
ni_N_link_d = 5  # the number of links directly downstream of this node
ni_Mlink_u1 = 6  # the index of the 1st link directly upstream
ni_Mlink_u2 = 7  # the index of the 2nd link directly upstream
ni_Mlink_u3 = 8  # the index of the 3rd link directly upstream
ni_Mlink_d1 = 9  # the index of the 1st link directly downstream
ni_Mlink_d2 = 10  # the index of the 2nd link directly downstream
ni_Mlink_d3 = 11  # the index of the 3rd link directly downstream
nr_directweight_u = 12  # the cumulative weight of the links directly upstream
nr_totalweight_u = 13  # the cumulative weight of all links upstream

l_left = 0  # number of links left in the list
link_id = 1  # link ID
li_idx = 2
li_link_type = 3  # int representation of the type of link
li_geometry = 4  # int representation of the type of link geometry
li_Mnode_u = 5  # ID of upstream node
li_Mnode_d = 6  # ID of downstream node
lr_Length = 7  # length of the link (0 if type != CONDUIT)
lr_Slope = 8  # average slope of the link, estimated with extreme points
lr_Roughness = 9  # Manning coefficient of the link (0 if type != CONDUIT)
lr_InitialFlowrate = 10  # initial flow rate
lr_InitialUpstreamDepth = 11  # initial upstream depth
lr_InitialDnstreamDepth = 12  # initial downstream depth

lr_Target = 1.0  # For the time being, the target length of an element is a hardcoded parameter

"""
Array Initialization:

This section defines the functions that initialize and populate the two governing arrays, nodes and links.

initialize_nodes_array():   we need to know the size of the array before we initialize it... but for some reason
                            file.seek(0) was not working with the csv.reader().  So I just wrote an independent
                            function that determines the size of the nodes_info.csv file and creates the nodes np.array
                            
populate_nodes_array(): with the nodes array initialized, this function reopens the nodes_info.csv file and, element
                        by element, populates the nodes array
                        
initialize_links_array(): same logic, but with a links array that is necessarily a different size than nodes

populate_links_array(): the corollary of populate_nodes_array, but for the links array

"""


def initialize_nodes_array(nodes_info):
    #print("Entering initialize_nodes_array()")
    global nodes, nrows
    with open(nodes_info, newline='') as csvfile:  # opens the file
        reader = csv.reader(csvfile, delimiter=',')  # defines the csv.reader, the delimiter is a comma
        nrows = -1
        for row in reader:  # count the number of rows in the file, we know it's square, so the number of columns is
                            # the length of any row
            nrows = nrows + 1
            ncols = len(row)
        multiprocessors = FileSettings.settingsdict['multiprocessors']
        nodes = np.empty((nrows + (multiprocessors - 1), ncols), dtype=object)  # excepting the the headers, initialize an np.empty array
                                                          # specify the dtype as object, which allows for any type
        """Adding empty rows of size multiprocessors to the end of the array ensures that there will be space for any 
        ghost nodes that are necessary.  For N multiprocessors, there can only be N ghost nodes.  Extraneous rows are
        removed at the end.  The -1 make sure the array (starting from 0) is the right size."""
    return


def populate_nodes_array(nodes_info):
    #print("Entering populate_nodes_array()")
    global nodes
    with open(nodes_info, newline='') as csvfile:  # opens the file
        reader = csv.reader(csvfile, delimiter=',')  # defines the csv.reader, the delimiter is a comma
        next(csvfile)  # omits the first line, which is the header information
        row_count = 0
        for row in reader:
            col_count = 0
            for col in row:
                nodes[row_count, col_count] = col  # for each (row, col) pair in the csv.reader, populate the
                                                   # corresponding location in nodes[row, col] with the value
                col_count = col_count + 1
            row_count = row_count + 1
    for i in range(len(nodes)):
        for j in range(len(nodes[1])):
            if nodes[i, j] is None:
                nodes[i, j] = -998877
    return


def initialize_links_array(links_info):
    #print("Entering initialize_links_array()")
    global links, lrows
    with open(links_info, newline='') as csvfile:  # opens the file
        reader = csv.reader(csvfile, delimiter=',')  # defines the csv.reader, the delimiter is a comma
        lrows = -1
        for row in reader:  # count the number of rows in the file, we know it's square, so the number of columns is
                            # the length of any row
            lrows = lrows + 1
            lcols = len(row)
        multiprocessors = FileSettings.settingsdict['multiprocessors']
        links = np.empty((lrows + (multiprocessors - 1), lcols), dtype=object)  # excepting the the headers, initialize an np.empty array
                                                          # specify the dtype as object, which allows for any type
    return


def populate_links_array(links_info):
    #print("Entering populate_links_array()")
    global links
    with open(links_info, newline='') as csvfile:  # opens the file
        reader = csv.reader(csvfile, delimiter=',')  # defines the csv.reader, the delimiter is a comma
        next(csvfile)  # omits the first line, which is the header information
        row_count = 0;
        for row in reader:
            col_count = 0
            for col in row:
                links[row_count, col_count] = col  # for each (row, col) pair in the csv.reader, populate the
                                                     # corresponding location in links[row, col] with the value
                col_count = col_count + 1
            row_count = row_count + 1
    for i in range(len(links)):
        for j in range(len(links[1])):
            if links[i, j] is None:
                links[i, j] = -998877
    return


def array_maker(nodes_info, links_info):
    initialize_nodes_array(nodes_info)
    populate_nodes_array(nodes_info)
    initialize_links_array(links_info)
    populate_links_array(links_info)
    #print("Exiting array_maker()")
    return


"""
Node Weighting:

This section of code defines the functions that govern the weight assignment to nodes

weighting_function(link_length):    the weight attributed to each link (that will ultimately be assigned to the
                                    downstream node) are normalized by lr_Target.  This gives an estimate of
                                    computational complexity. In the future lr_Target can be customized for each link.

null_value_convert(list):   this function is used to convert the null values that are defaulted in the nr_directweight_u
                            and nr_totalweight_u columns of the nodes array into float zeros.

local_node_weighting(): this function takes each node index (ni_idx) and finds that ni_idx in the downstream node column
                        of the links array (li_Mnode_d).  From this row in links the link weight is grabbed and ascribed
                        to the node-in-question's local weight (nr_directweight_u).

find_weight_index(root):    this function identifies the row that is being updated.  It is used as a row constant in the
                            upstream_weight_calculation function.

upstream_weight_calculation(root):  this is where the money is made.  This is a recursive method that maps through each
                                    root index in the nodes array, finds the upstream link index, jumps to that index in
                                    the link array, finds the upstream node for that link, and jumps to that.  All the 
                                    while the method is adding the local node weight to a running sum contained in
                                    nodes[weight_index, nr_totalweight_u]. This upstream traversing is called 
                                    recursively until a node is met that doesn't have an upstream link. The end result 
                                    of this function is that each node has been assigned essentially it's contributing 
                                    subsystem's 'complexity'.  The fact that this accumulating complexity is calculated
                                    without regard to the hydraulics of the systems gives the method incredible
                                    versatility.  The only prerequisite is an apriori definition of system direction.

nr_totalweight_assigner():  this function calls the upstream_weight_calculation() function for each node in the nodes
                            array.  It re-initializes the visited_flag list as all False before the recursive call, and
                            sets the new weight_index that will be receiving the weight updates within the function. It
                            also assigns the max_weight variable, which is later used as a factor in the threshold
                            decision variable.

find_endpoint():    This simple function just determines the number of nodes in a system that have no upstream neighbors
                    and therefore have a weight of zero.


**The two functions that need to be called in the Main module are local_node_weighting() to determine the direct weight
    applied to each node, and nr_totalweight_assigner() to determine the cumulative upstream weight for each node
"""


def weighting_function(link_length):  # The length of the link being weighted is passed from the links array
    global lr_Target                  # A constant for now, but a link attribute in the future
    return link_length/lr_Target      # The weight ascribed to the node is the number of expected elements on that link


def null_value_convert(list):             # Accepts a list argument
    #print("Entering null_value_convert()")
    for element in range(len(list)):      # Iterates through the list
        if list[element] == '-998877':    # If the element of the list is the default null value
            list[element] = 0.0           # Change that element to a float zero
    return


def local_node_weighting():
    #print("Entering local_node_weighting()")
    global lr_Target
    null_value_convert(nodes[:, nr_directweight_u])
    for row in range(len(nodes)):  # This weighting function occurs for each node
        rootnode_index = nodes[row, ni_idx]  # Assign to variable the node index
        links_row = 0  # Initialize the links_row which points to the index of the link
        for downstreamnode_index in links[:, li_Mnode_d]:  # Look in the downstream nodes column of the links array for
                                                            # the rootnode_index
            if downstreamnode_index == rootnode_index:  # When the rootnoode_index is found in that column
                nodes[row, nr_directweight_u] = float(nodes[row, nr_directweight_u]) + \
                    weighting_function(float(links[links_row, lr_Length]))
                                                                # Assign to the directweight column in the nodes array
                                                                # the sum of the weighting_function returned for each
                                                                # link directly upstream of that node
            links_row = links_row + 1  # Iterate the links_row, which ensures that the index in the link array being
                                       # weighed is consistent with the one that has the rootnode as a downstream node
    return


def find_weight_index(root_idx):    # Passes the root id that will set the row in the nodes array receiving the
                                                # running sum of upstream weights
    #print("Entering find_weight_index()")
    global weight_index
    for row in range(len(nodes)):               # For each row in the nodes array...
        root_identity = nodes[row, ni_idx]      # Assign to root_identity the node id from that row
        if root_idx == root_identity:           # If the root argument passed exists in this row of the nodes array...
            weight_index = row                  # Assign that row to a weight_index variable that will be static
                                                # throughout the recursive method assigning weights to that row
    return weight_index


def upstream_weight_calculation(root):
    """
    This function is the meat of the weighting calculation for each node.  Called recursively, this function goes
    through several loops.
    Loop 1: for each row in the nodes array
                if the nodes_id column matches the root AND that root has not been visited already on either the current
                call of upstream_weight_calculation OR other nr_totalweight_assigner() calls.
                **This last mask check ensures that weight of already partitioned subsystems do not contribute to the
                current network being assessed.
                    -change the visited_flag_weight for the root node to True so that it's weight cannot contribute more
                    than once to the weight_index
                    -update the nr_totalweight_u column in the weight_index (static) row of the nodes array with the
                    direct weight into the root node
                Loop 2: for each row in the links array
                        -look in the downstream endpoint of that link
                            if the root node index is found
                                -store the upstream endpoint node for that link
                            Loop 3: for each row in the nodes array once more
                                        if the upstream node is found
                                            -make that ni_idx the new root
                                            -recursively call the function on the new root
    visited_flag_weight --> boolean list to prevent double-counting links for a given subsystem traversal
    visited_network_mask --> boolean list to designate whether a node (and subsequently it's upstream subsystem)
                             has already been partitioned out of the system
    """
    global weight_index
    for node_row_number in range(len(nodes)):
        node_row_contents = nodes[node_row_number, ni_idx]
        if node_row_contents == -998877:
            continue
        if (root == node_row_contents) and (visited_flag_weight[node_row_number] is False) and \
                (visit_network_mask[node_row_number] is False):
            #root_idx = nodes[node_row_number, ni_idx]
            visited_flag_weight[node_row_number] = True
            #print("Weight index: " + str(weight_index) + " Node: " + str(nodes[weight_index, node_id]) + " Current Weight: " + str(nodes[weight_index, nr_totalweight_u]) + " Direct Weight: " + str(nodes[node_row_number, nr_directweight_u]))
            nodes[weight_index, nr_totalweight_u] = nodes[weight_index, nr_totalweight_u] + \
                                                    nodes[node_row_number, nr_directweight_u]
            #print("Row number: " + str(node_row_number) + ", Node ID: " + root + ", Node idx: " + root_idx)
            for link_row_number in range(len(links)):
                link_row_contents = links[link_row_number, li_Mnode_d]
                if node_row_contents == link_row_contents:
                    link_idx = links[link_row_number, li_idx]
                    node_upstream = links[link_row_number, li_Mnode_u]
                    #print("Node idx: " + root_idx + " is downstream of Link idx: " + link_idx + " which connects to Node idx: " + node_upstream)
                    for row_node_upstream in range(len(nodes)):
                        if node_upstream == nodes[row_node_upstream, ni_idx]:
                            new_root = node_upstream
                            upstream_weight_calculation(new_root)
    return


def nr_totalweight_assigner():
    print("Entering nr_totalweight_assigner()")
    """
    This function computes the cumulative upstream weight for each node in the system.  It does so by,
        -Initializing a new False mask/list for each new node being assigned weight.  This ensures that nodes that may
            exist on multiple paths upstream of a root node are not counted twice.
        -Determining the index of that primary node in the nodes array.  This index acts as the collector of the sum of
            upstream weights
        -Calling the recursive upstream_weight_calculation() function for that node index

    The max_weight variable is set as the maximum of the cumulative weights for each node
    """
    global visited_flag_weight, visit_network_mask, max_weight
    max_weight = 0
    for root in range(len(nodes)):
        if nodes[root, ni_idx] != -998877:
            visited_flag_weight = [False for node in nodes]
            weight_index = find_weight_index(nodes[root, ni_idx])
            upstream_weight_calculation(nodes[root, ni_idx])
    weight_list = [float(weight) for weight in nodes[:, nr_totalweight_u]]
    max_weight = max(weight_list)
    return


def find_endpoints():
    #print("Entering find_endpoints()")
    """
    This function is non-essential. Basically just determines if a node is an endpoint in the system by checking if it
    has any upstream links.
    """
    counter = 0
    mask_upstream_links = nodes[:, ni_Mlink_u1:ni_Mlink_u3+1]
    for row in mask_upstream_links:
        if all(up_link == '-998877' for up_link in row):
            counter = counter + 1
    return counter


"""
Subnetwork Identification and Optimal Partitioning:

Contains the functions that separate the nodes and lists arrays into multiple sub-arrays

Global Variables: subnetwork_container_nodes, subnetwork_container_links
    These variables are 3D zeros arrays that are initialized to be (depth x height x width) where 'depth' is the number
    of processors (i.e. the number of subsystems) to which the program will ascribe data, 'height' is the number of
    nodes or links, and 'width' contains the node/link info the same as the nodes, links arrays from the Array
    Initialization section.

subnetwork_carving(root):   given a root node_id, this function identifies the subnetwork that is upstream of this node
                            and copies that node's row information from the nodes array into subnetwork_container_nodes

subnetwork_links(): the subnetwork delineation is conducted on a node-wise basis of weights, but the vast majority of
                    the future elements will be contained within the links.  This function adds all of the links that
                    are "contained" within the subsystem to subnetworks_container_links.  Here, "contained" is defined
                    as having both endpoint nodes within the subsystem.

ideal_partition_check():    Computes the partition threshold and checks the network to determine if an optimal partition
                            is possible.  An optimal partition occurs if the partition threshold is found at its exact
                            value at a node.

All above functions can be used in the event that an optimal partition is possible for the network. Should an optimal
partition exist, the subnetwork_carving() routine is called to delineate the subnetwork upstream of that node.
Subsequently, the new network is re-weighted using local_node_weighting() and nr_totalweight_assigner().  Then the
ideal_partition_check is called to determine if optimal partitioning is still possible.  Each is called in main().
"""


def subnetwork_carving(root, i):
    """
    This function takes a node_id as an argument and determines the upstream subnetwork.  The traversal through the
    network is similar to the upstream_weight_calculation() function.

    For a node to be included in the root node's subnetwork, it must be found through the
    node_idx -> ni_Mlink_u -> li_idx -> li_Mnode_u -> node_idx
    pathway.

    If the visit_network_mask is False, it means that node has never been visited, and recursion can be called on
    that node. If the visit_network_mask is True, then that node already belongs to another subnetwork.  Furthermore,
    the entire upstream subnetwork of that node must already be accounted for.  To maintain connection between the
    subnetworks, the nodes with visit_network_mask == True are included in the current subnetwork, but not
    passed recursively.
    """
    global visit_network_mask, subnetwork_container_nodes
    for node_row_number in range(len(nodes)):
        node_row_contents = nodes[node_row_number, ni_idx]
        if (root == node_row_contents) and (visit_network_mask[node_row_number] is False):
            #root_idx = nodes[node_row_number, ni_idx]
            visit_network_mask[node_row_number] = True
            subnetwork_container_nodes[i, node_row_number] = nodes[node_row_number]
            #print("Row number: " + str(node_row_number) + ", Node ID: " + root + ", Node idx: " + root_idx)
            for link_row_number in range(len(links)):
                link_row_contents = links[link_row_number, li_Mnode_d]
                if node_row_contents == link_row_contents:
                    node_upstream = links[link_row_number, li_Mnode_u]
                    #print("Node idx: " + root_idx + " is downstream of Link idx: " + links[link_row_number, li_idx] + " which connects to Node idx: " + node_upstream)
                    for row_node_upstream in range(len(nodes)):
                        if node_upstream == nodes[row_node_upstream, ni_idx]:
                            new_root = node_upstream
                            subnetwork_carving(new_root, i)
        # This elif statement occurs if a node has already been added to another subnetwork. In this case, the node is
        # added to the subnetwork to maintain connectivity, but not passed recursively.
        elif (root == node_row_contents) and (visit_network_mask[node_row_number] is True):
            subnetwork_container_nodes[i, node_row_number] = nodes[node_row_number]
    return


def ideal_partition_check(processors_counter):
    #print("Entering ideal_partition_check()")
    """
    This function determines whether an optimal partition is possible.

    First, the partition_threshold is computed as a simple ratio between the max_weight of the current network and the
    number of processors remaining to accept data.

    The nodes array is searched for an index containing exactly this threshold value, returning the node_id if it is
    found.  Otherwise, the routine searches for the node with the nearest overestimate of the weight threshold.  The
    nearest overestimate weight is initialized as the max_weight, as this is guaranteed to both exist and be larger than
    the partition threshold.  In the same loop as the search for the ideal threshold, if the node's cumulative upstream
    weight value is greater than the threshold but less than the current nearest overestimate, the nearest overestimate
    is overwritten.  Should no ideal partition node be found the function exits, passing the node_id of the nearest
    overestimate.

    The end result of this function is that the ideal_exists check is either True or False.  If it is True, that means
    an ideal partition exists and the function passes the id of that root node. If it is False, then an ideal check does
    not exist and the nearest_overestimate_node is passed.  The need for this nearest overestimate will be explained in
    the Non-Ideal Partitioning section.
    """
    global max_weight, partition_threshold, ideal_exists
    ideal_root, nearest_overestimate_node = '', ''
    ideal_exists = False
    nearest_overestimate = max_weight
    # processors_left = FileSettings.settingsdict['multiprocessors'] - processors_counter
    # partition_threshold = max_weight/processors_left
    print(partition_threshold)
    for nodes_row in range(len(nodes)):
        # This if statement first checks if that node is ideal
        if nodes[nodes_row, nr_totalweight_u] == partition_threshold:
            ideal_root = nodes[nodes_row, ni_idx]
            ideal_exists = True
            return ideal_root
        # Otherwise, the nodes whose weights are greater than the partition threshold are searched.  If the weight is
        # greater than the partition_threshold AND less than the current nearest_overestimate, the nearest_overestimate
        # is overwritten.
        if (float(nodes[nodes_row, nr_totalweight_u]) > partition_threshold) and \
                (nodes[nodes_row, nr_totalweight_u] <= nearest_overestimate):
            nearest_overestimate = nodes[nodes_row, nr_totalweight_u]
            nearest_overestimate_node = nodes[nodes_row, ni_idx]
    return nearest_overestimate_node


"""
Non-Ideal Partitioning:

In the likely case that optimal partitioning is not possible (i.e. the partition threshold cannot be found exactly
at any node in the network), another partitioning strategy must be pursued.  In this exception strategy, the node whose
weight is the nearest overestimate of the partition threshold is passed as the root to subnetwork_carving(). This
subnetwork is used as a dummy network within which to search for a suitable partition. As a check, spanning_check() is
called to determine if the partition threshold is spanned by any link in the dummy network.

The node weights (which should maintain their values from the larger network) are searched this time for the
nearest underestimate of the partition threshold. Calling subnetwork_carving() on this node yields the first pass at
filling the current subnetwork to the appropriate size. The partition threshold is reduced by the weight on the root
node, and the dummy network is re-assigned weight to reflect the clipped out subnetwork component.

Once again use ideal_partition_check and spanning_check() to determine if the new weight is found exactly or else
spanned by a link.  If not, re-iterate the process of selecting the nearest underestimate node weight,
subnetwork_carving(), and re-assigning weights to the leftover nodes.

The cascade of decreasing-in-size subnetwork_carving() components are added to the subnetwork being delineated until
either the ideal_partition_check() or spanning_check() criteria are met to exactly fill the appropriate size.

spanning_check():   if an ideal partition does not exist it might be possible that a link spans the partition threshold.
                    'Spanning' here means that the weight at the upstream node is less than the threshold, but the sum
                    of the upstream node and a downstream link are greater than the threshold.  This indicates that a
                    phantom node can be placed along that link such that the weight at the phantom node is exactly the
                    value of the partition threshold.

linear_interpolator():  this simple function will be used to determine location of the phantom node along a link

phantom_naming_convenction():   when BIPQuick is translated into Fortran, its not possible to have different data types
                                cohabit the same array. For this reason, the phantom nodes and links will simply be
                                named in numerical order, starting from a value at least one order of magnitude greater
                                than the largest index.

phantom_node_generator():   this function amends the global nodes and links arrays adding 1 phantom node and 1 phantom
                            link to subdivide the link that spans the partition threshold.  This function uses
                            linear_interpolator() to determine the length of the phantom link, populates the *_id, *_idx
                            and upstream/downstream neighbor columns, while assigning default values to the other cells.
"""


def spanning_check(partition_threshold):
    #print("Entering spanning_check()")
    """
    This function is given a weighted network and a partition threshold value that does not exist at any node.  The
    function's job is to determine whether any link 'spans' the partition threshold, by establishing a range of weights.
    The minimum of this range is the weight value that exists at the upstream node, and the maximum of the range is
    given by the upstream node weight plus the weight contributed by the link itself, and exists just before the
    downstream node.

    An empty weight_range array is initialized to hold the minimum and maximum weights spanned by each link.  The
    weight_range array is populated, in position 0, for each link by the cumulative weight of the upstream node.  This
    represents the minimum cumulative weight found on that link. Weight_range position 1 is calculated by adding the
    weight of the link itself (using the weighting_function()) to the minimum weight in position 0.

    Once properly populated, the weight_array is compared against the partition threshold argument to determine if any
    links span the threshold.  If so, then that link index is returned.  If not, the function returns an empty string.
    """
    global weight_range, exactly_spanned
    spanning_link = ''
    exactly_spanned = False
    weight_range = np.zeros((len(links), 2))
    for link_row in range(len(links)):
        upstream_node = links[link_row, li_Mnode_u]
        for node_row in range(len(nodes)):
            if nodes[node_row, ni_idx] == upstream_node:
                # This adds the upstream node weight to the minimum weight column of weight_range
                weight_range[link_row, 0] = nodes[node_row, nr_totalweight_u]
        # This computes the weight from the link itself. The max weight column is link_weight + the min column
        link_weight = weighting_function(float(links[link_row, lr_Length]))
        weight_range[link_row, 1] = weight_range[link_row, 0] + link_weight
    print(weight_range)
    # Check if the partition_threshold is between the min-max weight for any link
    for weight_row in range(len(weight_range)):
        if (weight_range[weight_row, 0] < partition_threshold < weight_range[weight_row, 1]) and \
                partition_boolean[weight_row] is False:
            spanning_link = links[weight_row, li_idx]
            partition_boolean[weight_row] = True
            # Return a link index
            return spanning_link
    # Otherwise return an empty string
    return spanning_link


def linear_interpolator(partition_threshold, spanning_link):
    #print("Entering linear_interpolator()")
    """
    This function is given a partition threshold and the spanning link.  It determines how far from the upstream
    node the phantom node must be placed.

    The key to this function is that the ratio of weights will be maintained in a conversion to length space. I.e.
    the weight_ratio is computed from the linear interpolation equation,

    (x - x1)/(x2 - x1) = weight_ratio   where x = partition_threshold, x1 = min weight, x2 = max weight

    Then, the length from the upstream node along the link at which the partition threshold will be found can be
    computed by multiplying the total length of the link by this weight_ratio.  This value is then returned.
    """
    global weight_range
    for link_row in range(len(links)):
        if links[link_row, li_idx] == spanning_link:
            total_length = float(links[link_row, lr_Length])
            start = weight_range[link_row, 0]
    weight_ratio = (partition_threshold - start)/weighting_function(total_length)
    length_from_start = weight_ratio*total_length
    print("The length from upstream node is: " + str(length_from_start))
    return length_from_start


def phantom_naming_convention():
    #print("Entering phantom_naming_convention()")
    """
    This function serves to solve the problem of what to call the newly created phantom nodes and links.  It determines
    the maximum number of digits that exists in the index column of the nodes and links array.  The phantom ni_idx and
    corresponding li_idx will begin at "1" starting at the next order of magnitude above the maximum.
    """
    global phantom_idx
    # Find the largest number in the index column of the nodes array
    node_indices = [int(index) for index in nodes[:, ni_idx]]
    max_node_idx = max(node_indices)

    # Find the largest number in the index column of the links array
    link_indices = [int(index) for index in links[:, li_idx]]
    max_link_idx = max(link_indices)

    # Determine the number of digits in the largest index number in either array. Set the starting value of the phantom
    # nodes/links to be the next order of magnitude higher.
    digits = len(str(max(max_node_idx, max_link_idx)))
    phantom_idx = 10**digits
    return


def phantom_node_generator(spanning_link, length_from_start, i):
    #print("Entering phantom_node_generator()")
    """
    The purpose of this function is to amend the nodes and links arrays with a new node and link.  It is critical to
    note that this new node and link do not change the geometry of the system, only the topology. The new node is
    positioned such that the cumulative weight of the node is exactly the partition threshold.  The new link inherits
    most of the properties of the link it once belonged to, with the exception of the length, which has been split to
    reflect that one link has been separated into two.
    """
    global counter, partition_threshold, nodes, links, subnetwork_container_nodes, subnetwork_container_links, nrows, lrows, phantom_array_loc

    # The name of the phantom node is set as a string-itized integer, the location is known as the first "empty" row
    phantom_name = str(phantom_idx + counter)
    phantom_array_loc = nrows + i
    phantom_array_loc_link = lrows + i


    # For the last row of the new nodes array, all of the column values can be deduced.
    # This column doesn't matter at all, it just needs a value
    nodes[phantom_array_loc, nleft] = '0'
    # The phantom node id/idx is named to reflect how many phantom nodes exist with the integer counter suffix
    nodes[phantom_array_loc, node_id] = phantom_name
    nodes[phantom_array_loc, ni_idx] = phantom_name
    # Node Type '0' denotes a simple 1 in - 1 out junction, which all phantom nodes are guaranteed to be
    nodes[phantom_array_loc, ni_node_type] = '0'
    # The phantom nodes are also guaranteed to only have 1 upstream and 1 downstream link
    nodes[phantom_array_loc, ni_N_link_u] = '1'
    nodes[phantom_array_loc, ni_N_link_d] = '1'
    # Initialize the first upstream link map as the old spanning link, and the downstream link map as the phantom link
    nodes[phantom_array_loc, ni_Mlink_u2:ni_Mlink_d3+1] = '-998877'
    nodes[phantom_array_loc, ni_Mlink_u1] = spanning_link
    nodes[phantom_array_loc, ni_Mlink_d1] = phantom_name
    # The direct weight doesn't matter as it is the locus for a subnetwork partition
    nodes[phantom_array_loc, nr_directweight_u] = 0
    # The cumulative weight at the phantom node is also known.  By design it is the partition_threshold.
    nodes[phantom_array_loc, nr_totalweight_u] = partition_threshold

    # If this partition is not the last one, go ahead and add that node to the next partition also
    # if i != FileSettings.settingsdict['multiprocessors'] - 1:
    #     subnetwork_container_nodes[i+1, phantom_array_loc] = nodes[phantom_array_loc]

    # This code block is for updating the links array
    for link_row in range(len(links)):
        #new_links[link_row] = links[link_row]
        # Instead of adding an empty row, the phantom link is copied from the spanning link and modified
        if links[link_row, li_idx] == spanning_link:
            links[phantom_array_loc_link] = links[link_row]
            # The new length of the phantom link is the original length minus the length_from_start node location
            links[phantom_array_loc_link, lr_Length] = float(links[link_row, lr_Length]) - length_from_start
            # When the phantom node is inserted, the direct weight into the downstream node of the spanning link is
            # reduced.  The exact new value of the direct weight is computed using the weighting_function() on the
            # length of the induced phantom link
            downstream_node = links[link_row, li_Mnode_d]
            for node_row in range(len(nodes)):
                if nodes[node_row, ni_idx] == downstream_node:
                    nodes[node_row, nr_directweight_u] = nodes[node_row, nr_directweight_u] - weighting_function(length_from_start)
            # The downstream node map for the spanning link is updated to the phantom node and the length is updated
            links[link_row, li_Mnode_d] = str(phantom_idx + counter)
            links[link_row, lr_Length] = length_from_start
            break

    # The geometric properties of the phantom link match those from the spanning link.  As previously stated, the
    # geometry of the system is conserved, only the topology changes.
    links[phantom_array_loc_link, l_left] = '0'
    # The link id/idx is set similarly to the phantom node convention
    links[phantom_array_loc_link, link_id] = str(phantom_idx + counter)
    links[phantom_array_loc_link, li_idx] = str(phantom_idx + counter)
    #The upstream node map for the phantom link is set to be the phantom node
    links[phantom_array_loc_link, li_Mnode_u] = nodes[phantom_array_loc, ni_idx]
    #links = new_links

    # The next time a phantom node is required, it will be dubbed '102'
    counter = counter + 1

    return


"""
Post-Processing Module

The Main module ends with a populated 3D array of the nodes that are contained within each partition.  The Post-
Processing module much first match these node-boundary definitions with the links that reside between them.

This section of the code then transforms the 3D arrays containing the nodes and links information for each partition
back into a 2D array that contains the whole system (including phantom nodes/links but excluding placeholders).  This 2D
array will be arranged such that in an automatic n-ways splitting of the array will keep the elements on one partition
together.

subnetworks_links(nodes_container):  This function populates the subnetwork_container_links array with the links
                                        "contained" within a subnetwork.  A set of endpoints (upstream and downstream
                                        nodes; each link has exactly one of each) is created for each link.  If this
                                        set is a subset of the potential_endpoints that exist within the ni_idx column
                                        of the nodes_container, then both endpoints for a link are in the subnetwork.
                                        Therefore, the link must be contained.

reorganize_arrays(nodes_container, links_container): This function consolidates the 3D nodes/links containers into 2D
                                                     arrays are arranged in terms of the partition they appear in.  If
                                                     these arrays were to be passed to Fortran CoArrays, they would be
                                                     sent to processors in alignment with the partitioning scheme.
"""


def subnetworks_links(nodes_container):
    """
    This function works by iterating through the links array, and comparing the endpoints of each link to the
    nodes_container to determine if that link is completely contained within the partitioned subsystem.  If so, that
    link is added to the accounted_for_links array.  This array is maintained in the case of cross-connections in the
    system that are separated by a partition.  It could be that both nodes of a link appear in both partitions, but only
    the preceding partition is permitted access to that link.

    The links for each partition are added to a container in the same fashion as the nodes.
    """
    global subnetwork_container_links, i, accounted_for_links
    for link_row in range(len(links)):
        endpoints = links[link_row, li_Mnode_u:li_Mnode_d+1]
        potential_endpoints = nodes_container[:, ni_idx]
        # Set A is less than or equal to set be if all of the elements in set A are also contained within set B.
        # Need to check the previous nodes_container, because it might already have that link
        accounted_for_links = np.empty((len(links), 1), dtype=object)
        for past_partitions in range(i):
            link_counter = 0
            for link_row_contents in subnetwork_container_links[past_partitions, :, li_idx]:
                if link_row_contents is not None:
                    accounted_for_links[link_counter, 0] = link_row_contents
                link_counter = link_counter + 1
        if (set(endpoints) <= set(potential_endpoints)) and (links[link_row, li_idx] not in accounted_for_links):
            # For our case, the set is two elements.  If both elements are in the subsystem, then that link is contained
            subnetwork_container_links[i, link_row] = links[link_row]
    return


def populate_last_processor():
    """Add the nodes belonging to the links that haven't yet been accounted for to the final processor"""
    global subnetwork_container_nodes, subnetwork_container_links, nodes, links, accounted_for_links
    for link_row in range(len(links)):
        if links[link_row, li_idx] not in accounted_for_links:
            subnetwork_container_links[-1, link_row] = links[link_row]
    for link_row in range(len(subnetwork_container_links[-1])):
        for node_row in range(len(nodes)):
            if subnetwork_container_links[-1, link_row, li_Mnode_u] == nodes[node_row, ni_idx]:
                subnetwork_container_nodes[-1, node_row] = nodes[node_row]
            if subnetwork_container_links[-1, link_row, li_Mnode_d] == nodes[node_row, ni_idx]:
                subnetwork_container_nodes[-1, node_row] = nodes[node_row]
    return


def reorganize_arrays(nodes_container, links_container):
    """
    The purpose of this function is to take the n-processor length containers for the nodes and links that appear on
    each partition, and recombine them each into a single array that is the ultimate return of the BIPQuick routine.
    """
    global nodes, links
    # Initialize two 2D arrays that are the same size as the nodes/links arrays
    reorganized_nodes = np.empty((len(nodes), len(nodes[1])), dtype=object)
    reorganized_links = np.empty((len(links), len(links[1])), dtype=object)
    # Initialize node/link counters which will be used to ensure that the
    nodes_row_counter = 0
    links_row_counter = 0
    for i in range(len(nodes_container)):
        for row in nodes_container[i]:
            if None in row:
                continue
            elif row[ni_idx] in reorganized_nodes[:, ni_idx]:
                continue
            elif row[ni_idx] == -998877:
                continue
            else:
                reorganized_nodes[nodes_row_counter] = row
                nodes_row_counter = nodes_row_counter + 1
    clipped_nodes = np.empty((nodes_row_counter, len(nodes[1])), dtype=object)
    for row in range(len(clipped_nodes)):
        clipped_nodes[row] = reorganized_nodes[row]
    nodes = clipped_nodes

    for i in range(len(links_container)):
        for row in links_container[i]:
            if None in row:
                continue
            elif row[ni_idx] in reorganized_links[:, ni_idx]:
                continue
            elif row[ni_idx] == -998877:
                continue
            else:
                reorganized_links[links_row_counter] = row
                links_row_counter = links_row_counter + 1
    clipped_links = np.empty((links_row_counter, len(links[1])), dtype=object)
    for row in range(len(clipped_links)):
        clipped_links[row] = reorganized_links[row]
    links = clipped_links
    return


def post_processing():
    global subnetwork_container_nodes, subnetwork_container_links, i
    for i in range(len(subnetwork_container_nodes)):
        nodes_container = subnetwork_container_nodes[i]
        subnetworks_links(nodes_container)
        #populate_last_processor()
    reorganize_arrays(subnetwork_container_nodes, subnetwork_container_links)
    return



main()
print(subnetwork_container_nodes)
post_processing()
print(subnetwork_container_links)
print(nodes)
print(links)

