import logging
from Bio.Seq import Seq

class Feature():
    SO_REGION = 'http://identifiers.org/so/SO:0000001'
    SO_SEQUENCE_FEATURE = 'http://identifiers.org/so/SO:0000110'

    GENERIC_ROLES = {
        SO_REGION,
        SO_SEQUENCE_FEATURE
    }

    def __init__(self, nucleotides, identity, roles, sub_identities=[], parent_identities=[]):
        self.nucleotides = nucleotides
        self.identity = identity
        self.sub_identities = sub_identities
        self.parent_identities = parent_identities
        self.roles = set(roles)

        self.logger = logging.getLogger('synbict')

    def reverse_complement_nucleotides(self):
        return str(Seq(self.nucleotides).reverse_complement())

    @classmethod
    def has_non_generic_role(cls, roles):
        return len(roles.difference(cls.GENERIC_ROLES)) > 0

    def is_non_generic(self):
        return self.has_non_generic_role(self.roles)