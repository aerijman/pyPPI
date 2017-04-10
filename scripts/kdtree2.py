#!/usr/bin/python
# encoding: utf-8

""" KDTree implementation.

Features:

- nearest neighbours search

Matej Drame [matej.drame@gmail.com]
"""

__version__ = "1r11.1.2010"
__all__ = ["KDTree"]
import heapq

def square_distance(pointA, pointB):
    # squared euclidean distance
    distance = 0.0
    dimensions = len(pointA) # assumes both points have the same dimensions
    for dimension in range(dimensions):
        distance += (pointA[dimension] - pointB[dimension])**2
    return distance

class KDTreeNode():
    def __init__(self, node, left, right):
        self.point = node.coord
        self.node = node
        self.left = left
        self.right = right
    
    def is_leaf(self):
        return (self.left == None and self.right == None)
       
class KDTree():
    """ KDTree implementation.
    
        Example usage:
        
            from kdtree import KDTree
            
            data = <load data> # iterable of points (which are also iterable, same length)
            point = <the point of which neighbours we're looking for>
            
            tree = KDTree.construct_from_data(data)
            nearest = tree.query(point, t=4) # find nearest 4 points
    """
    
    def __init__(self, data):
        def build_kdtree(point_list, depth):
            # code based on wikipedia article: http://en.wikipedia.org/wiki/Kd-tree
            if not point_list:
                return None

            # select axis based on depth so that axis cycles through all valid values
            axis = depth % len(point_list[0].coord) # assumes all points have the same dimension

            # sort point list and choose median as pivot point,
            # TODO: better selection method, linear-time selection, distribution
            point_list.sort(key=lambda node: node.coord[axis])
            median = len(point_list)/2 # choose median

            # create node and recursively construct subtrees
            node = KDTreeNode(node=point_list[median],
                              left=build_kdtree(point_list[0:median], depth+1),
                              right=build_kdtree(point_list[median+1:], depth+1))
            return node
        
        self.root_node = build_kdtree(data, depth=0)
    
    @staticmethod
    def construct_from_data(data):
        tree = KDTree(data)
        return tree

    """
    distance - pow distance
    
    """
    def findByDistance(self, query_point,distance):
        if self.root_node == None:
            return
        
        #statistics = {'nodes_visited': 0, 'far_search': 0, 'leafs_reached': 0}
        res = []
        def conditionalAdd(node):
            dist=square_distance(node.point,query_point)
            if dist<distance:
                return True
                res.append(node.node)

        def nn_search(node, query_point, distance, depth, lastDist=0):
            if node == None:
                return
            
            #statistics['nodes_visited'] += 1
            
            
            if node.is_leaf():#and square_distnace moved below... 21/5/13
                #statistics['leafs_reached'] += 1
                if square_distance(node.point,query_point)<distance:
                    yield node.node
                return
            
            axis = depth % len(query_point)
            
            # figure out which subtree to search
            near_subtree = None # near subtree
            far_subtree = None # far subtree (perhaps we'll have to traverse it as well)
            
            # compare query_point and point of current node in selected dimension
            # and figure out which subtree is farther than the other
            if query_point[axis] < node.point[axis]:
                near_subtree = node.left
                far_subtree = node.right
            else:
                near_subtree = node.right
                far_subtree = node.left
            
            for r in nn_search(near_subtree, query_point, distance, depth+1):
                yield r
            currDist=(node.point[axis] - query_point[axis])**2
            if currDist+lastDist< distance:
                if square_distance(node.point,query_point)<distance:
                    yield node.node
                #conditionalAdd(node)
                #statistics['far_search'] += 1
                for r in nn_search(far_subtree, query_point, distance, depth+1,lastDist=currDist):
                    yield r
            return

        
        
        for r in nn_search(self.root_node, query_point, distance, depth=0):
            yield r

    """
    Finds nearest atoms to query point
    """
    def findNearest(self, query_point, num):
        if self.root_node == None:
            return
        res = []
        def conditionalAdd(node):
            distance = square_distance(node.point,query_point)
            if distance == 0.0:
                distance = 0.000000001
            if len(res)==num:
                heapq.heappushpop(res, (1.0/distance, distance, node.node))
            else:
                heapq.heappush(res, (1.0/distance, distance, node.node))
            
        def nn_search(node, query_point, depth, lastDist=0):
            if node == None:
                return
            
            if node.is_leaf():
                #statistics['leafs_reached'] += 1
                conditionalAdd(node)
                return
            
            axis = depth % len(query_point)
            
            # figure out which subtree to search
            near_subtree = None # near subtree
            far_subtree = None # far subtree (perhaps we'll have to traverse it as well)
            
            # compare query_point and point of current node in selected dimension
            # and figure out which subtree is farther than the other
            if query_point[axis] < node.point[axis]:
                near_subtree = node.left
                far_subtree = node.right
            else:
                near_subtree = node.right
                far_subtree = node.left

            
            nn_search(near_subtree, query_point, depth+1)

            currDist=(node.point[axis] - query_point[axis])**2
            if len(res)<num or currDist+lastDist< res[0][1]:
                conditionalAdd(node)
                nn_search(far_subtree, query_point, depth+1,lastDist=currDist)
            return

        nn_search(self.root_node, query_point, depth=0)

        sorted(res,key=lambda a:a[1])
        return [x[2] for x in res]

    def findFarest(self, query_point, num):
        if self.root_node == None:
            return

        res=[]
        def conditionalAdd(node):
            distance=square_distance(node.point,query_point)
            if len(res)==num:
                heapq.heappushpop(res, (distance, node.node))
            else:
                heapq.heappush(res, (distance, node.node))
                
        def nn_search(node, query_point, depth, lastDist=0):
            
            if node == None:
                return
            
            
            if node.is_leaf():
                #statistics['leafs_reached'] += 1
                conditionalAdd(node)
                return
            
            axis = depth % len(query_point)
            
            # figure out which subtree to search
            near_subtree = None # near subtree
            far_subtree = None # far subtree (perhaps we'll have to traverse it as well)
            
            # compare query_point and point of current node in selected dimension
            # and figure out which subtree is farther than the other
            if query_point[axis] < node.point[axis]:
                near_subtree = node.left
                far_subtree = node.right
            else:
                near_subtree = node.right
                far_subtree = node.left
            
            nn_search(far_subtree, query_point, depth+1)
            currDist=(node.point[axis] - query_point[axis])**2
            if len(res)<num or currDist+lastDist> res[0][1]:
                conditionalAdd(node)
                nn_search(near_subtree, query_point, depth+1,lastDist=currDist)
                    
            
        
        nn_search(self.root_node, query_point, depth=0)
        sorted(res,key=lambda a:-a[0])
        return [x[1] for x in res]

