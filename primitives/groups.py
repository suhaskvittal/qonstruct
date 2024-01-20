"""
    author: Suhas Vittal
    date:   17 January 2024
"""

import math

class GroupElement:
    """
        Implementation of a group element. Two properties:
            (1) name (how to represent the element)
            (2) value (i.e. an integer).
        mul_operator can also be set to define multiplication between group elements.

        If mul_operator=None, then it is assumed that "value" defines *.
        If set, mul_operator should be a lambda that takes in two elements and returns one.
    """
    def __init__(self, name: str, value: any, mul_operator=None):
        self.name = name
        self.value = value
        # Set muliplication operator.
        if mul_operator is None:
            self.mul_operator = __mul_default__
        else:
            self.mul_operator = mul_operator

    def __eq__(self, other: GroupElement) -> bool:
        return self.value == other.value
    
    def __hash__(self) -> int:
        return hash(self.value)
    
    def __mul__(self, other: GroupElement) -> GroupElement:
        return mul_operator(self, other)

    def __pow__(self, p: int) -> GroupElement:
        return mul_operator(self, self)

    def __repr__(self) -> str:
        return self.name

    def __mul_default__(a: GroupElement, b: GroupElement):
        prod = GroupElement('%s.%s' % (a.name, b.name), a*b)
        return prod

class Group:
    def __init__(self, identity: GroupElement, generators: list[GroupElement]):
        self.identity = identity;
        self.generators = generators
        # Memoizer data structures:
        self._order_map = {}
        self._inv_map = {}
        # Now, go ahead and compute the group order.
        self.order = math.lcm(order(g) for g in generators)

    def contains(

    def inv(self, x: GroupElement) -> GroupElement:
        if x == identity:
            return identity
        else:
            return x ** (order(x)-1)

    def order(self, x: GroupElement) -> int:
        if x in self._order_map:
            return x
        # Else, we need to compute the order if it is not memoized.
        LIMIT_ASSUME_ORDER_INFINITE = 1_000_000
        curr = x
        r = 1
        while curr != identity && r < LIMIT_ASSUME_ORDER_INFINITE:
            curr = curr * x
            r += 1
        self._order_map[x] = r
        return r

