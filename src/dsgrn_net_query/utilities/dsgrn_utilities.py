def satisfies_hex_constraints(param, hex_constraints):
    # param is a DSGRN.Parameter object
    #
    # hex_constraints is a dictionary of lists keying a tuple of two integers to a list
    # of hex numbers.
    # The tuple key describes the node type: (num inedges, num outedges) and the list contains the
    # allowable hex codes for this node type. In the algorithm below, only those hex codes in the list
    # are permitted for the node type.
    # Example: {(1,2) : ["C", "8"], (3,1) : ["0"]} means that any node with 1 in-edge and 2 out-edges must have
    # hex code 0x0C or 8 and any node with 3 in-edges and 1 outedge must have hex code 0.
    #
    # Node types not listed in hex_constraints are ignored.
    # Node types with different algebraic forms are not distinguished. For example, a node with input A+B+C and one
    # out-edge and a node with (A+B)(~C) and one out-edge are both called (3,1).
    # TODO: Extend API to distinguish between logic and type of regulation (activating or repressing).

    node_types = [(v, (v.numInputs(), v.numOutputs())) for v in param.logic()]
    for v,n in node_types:
        if n in hex_constraints and v.hex() not in hex_constraints[n]:
            return False
    return True










