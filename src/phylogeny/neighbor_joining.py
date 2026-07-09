"""Reusable Neighbor Joining tree reconstruction from a labelled distance matrix."""

from __future__ import annotations

from phylogeny.upgma import Cluster, validate_distance_matrix


DistanceMap = dict[tuple[str, str], float]


def neighbor_joining(labels: list[str], matrix: list[list[float]]) -> Cluster:
    """Reconstruct a distance-based tree with the Neighbor Joining algorithm.

    :param labels: Ordered taxon labels corresponding to matrix rows and columns.
    :param matrix: Square symmetric matrix of pairwise taxon distances.
    :return: Root cluster of the reconstructed tree.
    """
    validate_distance_matrix(labels, matrix)

    # Keep active clusters sorted by deterministic names so Q-matrix ties resolve reproducibly.
    clusters = [
        Cluster(name=label, members=(label,), height=0.0)
        for label in sorted(labels)
    ]
    distances = _initial_distance_map(labels, matrix)

    while len(clusters) > 3:
        q_values = q_matrix(clusters, distances)
        first, second = select_neighbor_pair(q_values, clusters)
        left, right = sorted((first, second), key=lambda cluster: cluster.name)
        branch_left, branch_right = neighbor_branch_lengths(left, right, clusters, distances)

        # Internal names are built from sorted descendant labels to keep Newick output stable.
        merged_members = tuple(sorted(left.members + right.members))
        merged = Cluster(
            name="+".join(merged_members),
            members=merged_members,
            height=0.0,
            left=left,
            right=right,
            left_branch_length=branch_left,
            right_branch_length=branch_right,
        )

        distances = updated_distance_map(left, right, merged, clusters, distances)
        clusters = [
            cluster for cluster in clusters
            if cluster is not first and cluster is not second
        ]
        clusters.append(merged)
        clusters.sort(key=lambda cluster: cluster.name)

    if len(clusters) == 2:
        # A two-taxon input has no NJ Q-matrix step; connect both leaves at the midpoint.
        left, right = sorted(clusters, key=lambda cluster: cluster.name)
        final_distance = lookup_distance(left, right, distances)
        return Cluster(
            name="+".join(sorted(left.members + right.members)),
            members=tuple(sorted(left.members + right.members)),
            height=0.0,
            left=left,
            right=right,
            left_branch_length=final_distance / 2.0,
            right_branch_length=final_distance / 2.0,
        )

    # Standard NJ finishes with three clusters attached to one unrooted central node.
    left, right, third = sorted(clusters, key=lambda cluster: cluster.name)
    left_length, right_length, third_length = final_branch_lengths(
        left,
        right,
        third,
        distances,
    )
    return Cluster(
        name="+".join(sorted(left.members + right.members + third.members)),
        members=tuple(sorted(left.members + right.members + third.members)),
        height=0.0,
        left=left,
        right=right,
        third=third,
        left_branch_length=left_length,
        right_branch_length=right_length,
        third_branch_length=third_length,
    )


def q_matrix(clusters: list[Cluster], distances: DistanceMap) -> DistanceMap:
    """Construct the Neighbor Joining Q matrix for the active clusters.

    :param clusters: Current active clusters in the reconstruction.
    :param distances: Pairwise distances between active clusters.
    :return: Mapping of unordered cluster name pairs to their Q-matrix value.
    """
    if len(clusters) < 3:
        raise ValueError("Q matrix requires at least three active clusters")

    # Totals are reused for every Q entry in the current iteration.
    totals = total_distances(clusters, distances)
    active_count = len(clusters)
    return {
        _distance_key(left.name, right.name): (active_count - 2)
        * lookup_distance(left, right, distances)
        - totals[left.name]
        - totals[right.name]
        for left_index, left in enumerate(clusters)
        for right in clusters[left_index + 1:]
    }


def total_distances(clusters: list[Cluster], distances: DistanceMap) -> dict[str, float]:
    """Calculate each active cluster's total distance to every other cluster.

    :param clusters: Current active clusters in the reconstruction.
    :param distances: Pairwise distances between active clusters.
    :return: Mapping from cluster name to summed distance.
    """
    totals: dict[str, float] = {}
    for cluster in clusters:
        # Each total excludes the cluster itself and sums all current distances from it.
        totals[cluster.name] = sum(
            lookup_distance(cluster, other, distances)
            for other in clusters
            if other is not cluster
        )
    return totals


def select_neighbor_pair(
    q_values: DistanceMap,
    clusters: list[Cluster],
) -> tuple[Cluster, Cluster]:
    """Select the pair with the smallest Q value using deterministic tie-breaking.

    :param q_values: Q-matrix values keyed by unordered cluster name pair.
    :param clusters: Current active clusters in the reconstruction.
    :return: Selected pair of neighboring clusters.
    """
    if not q_values:
        raise ValueError("Q matrix must contain at least one candidate pair")

    # Tied Q values use alphabetical cluster names, matching the UPGMA tie-breaking style.
    selected_names, _ = min(
        q_values.items(),
        key=lambda item: (item[1], item[0]),
    )
    clusters_by_name = {cluster.name: cluster for cluster in clusters}
    try:
        return clusters_by_name[selected_names[0]], clusters_by_name[selected_names[1]]
    except KeyError as error:
        raise ValueError(f"Q matrix contains unknown cluster: {error.args[0]}") from error


def neighbor_branch_lengths(
    first: Cluster,
    second: Cluster,
    clusters: list[Cluster],
    distances: DistanceMap,
) -> tuple[float, float]:
    """Calculate branch lengths from a new NJ node to its two joined children.

    :param first: First child cluster in deterministic left/right order.
    :param second: Second child cluster in deterministic left/right order.
    :param clusters: Current active clusters before joining the pair.
    :param distances: Pairwise distances between active clusters.
    :return: Branch lengths to first and second, respectively.
    """
    if len(clusters) < 3:
        raise ValueError("Neighbor Joining branch lengths require at least three clusters")

    totals = total_distances(clusters, distances)
    pair_distance = lookup_distance(first, second, distances)
    divisor = 2.0 * (len(clusters) - 2)
    # NJ allows child branches to differ when one lineage is farther from the remaining taxa.
    first_length = 0.5 * pair_distance + (totals[first.name] - totals[second.name]) / divisor
    second_length = pair_distance - first_length
    return first_length, second_length


def updated_distance_map(
    first: Cluster,
    second: Cluster,
    merged: Cluster,
    clusters: list[Cluster],
    distances: DistanceMap,
) -> DistanceMap:
    """Update active cluster distances after joining one NJ pair.

    :param first: First joined child cluster.
    :param second: Second joined child cluster.
    :param merged: Newly created parent cluster.
    :param clusters: Active clusters before removing the joined children.
    :param distances: Pairwise distances between active clusters.
    :return: Distance map for the next NJ iteration.
    """
    next_distances: DistanceMap = {}
    remaining = [
        cluster for cluster in clusters
        if cluster is not first and cluster is not second
    ]

    for left_index, left in enumerate(remaining):
        for right in remaining[left_index + 1:]:
            # Distances between unaffected clusters carry forward unchanged.
            next_distances[_distance_key(left.name, right.name)] = lookup_distance(
                left,
                right,
                distances,
            )

    for cluster in remaining:
        # NJ updates distances without cluster-size weighting.
        updated = (
            lookup_distance(first, cluster, distances)
            + lookup_distance(second, cluster, distances)
            - lookup_distance(first, second, distances)
        ) / 2.0
        next_distances[_distance_key(merged.name, cluster.name)] = updated

    return next_distances


def final_branch_lengths(
    first: Cluster,
    second: Cluster,
    third: Cluster,
    distances: DistanceMap,
) -> tuple[float, float, float]:
    """Calculate branch lengths for the final three-way NJ node.

    :param first: First final cluster in deterministic output order.
    :param second: Second final cluster in deterministic output order.
    :param third: Third final cluster in deterministic output order.
    :param distances: Pairwise distances between the final active clusters.
    :return: Branch lengths to first, second and third, respectively.
    """
    first_second = lookup_distance(first, second, distances)
    first_third = lookup_distance(first, third, distances)
    second_third = lookup_distance(second, third, distances)

    # With three clusters left, each branch is determined by the three pairwise distances.
    first_length = (first_second + first_third - second_third) / 2.0
    second_length = (first_second + second_third - first_third) / 2.0
    third_length = (first_third + second_third - first_second) / 2.0
    return first_length, second_length, third_length


def lookup_distance(first: Cluster, second: Cluster, distances: DistanceMap) -> float:
    """Return the stored distance between two active clusters.

    :param first: First cluster in the pair.
    :param second: Second cluster in the pair.
    :param distances: Pairwise distances keyed by cluster name pair.
    :return: Stored pairwise distance.
    """
    try:
        # Distances are stored once per unordered pair.
        return distances[_distance_key(first.name, second.name)]
    except KeyError as error:
        raise ValueError(
            f"Missing distance between clusters {first.name!r} and {second.name!r}"
        ) from error


def _initial_distance_map(labels: list[str], matrix: list[list[float]]) -> DistanceMap:
    """Build the active distance map from the validated input matrix.

    :param labels: Ordered taxon labels corresponding to matrix rows and columns.
    :param matrix: Square symmetric matrix of pairwise taxon distances.
    :return: Pairwise distances keyed by unordered taxon label pairs.
    """
    label_indexes = {label: index for index, label in enumerate(labels)}
    sorted_labels = sorted(labels)
    return {
        _distance_key(left, right): matrix[label_indexes[left]][label_indexes[right]]
        for left_index, left in enumerate(sorted_labels)
        for right in sorted_labels[left_index + 1:]
    }


def _distance_key(first: str, second: str) -> tuple[str, str]:
    """Return a stable key for an unordered cluster-name pair.

    :param first: First cluster name.
    :param second: Second cluster name.
    :return: Alphabetically ordered pair of names.
    """
    if first == second:
        raise ValueError("A distance key requires two distinct cluster names")
    return tuple(sorted((first, second)))
