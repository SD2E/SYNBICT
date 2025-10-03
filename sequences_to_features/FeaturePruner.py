import sbol2
import Feature
import logging

class FeaturePruner():

    COMMON_ROLE_DICT = {
        sbol2.SO_PROMOTER: 'promoter',
        sbol2.SO_CDS: 'CDS',
        sbol2.SO_TERMINATOR: 'terminator',
        'http://identifiers.org/so/SO:0001977': 'ribonuclease_site'
    }

    def __init__(self, feature_library, roles=set()):
        self.feature_library = feature_library
        self.roles = roles

        self.logger = logging.getLogger('synbict')

    @classmethod
    def __has_min_length(cls, feature, min_feature_length):
        return min_feature_length == 0 or len(feature.nucleotides) >= min_feature_length

    @classmethod
    def __is_covered(cls, anno, cover_annos, cover_offset):
        for cover_anno in cover_annos:
            if not abs(cover_anno[0] - anno[0]) <= cover_offset or not abs(cover_anno[1] - anno[1]) <= cover_offset:
                return False

        return True

    def __remove_annotations(self, indices, annos, target_definition):
        for i in range(len(annos) - 1, -1, -1):
            if annos[i][5] is None:
                feature_identity = annos[i][2]
            else:
                feature_identity = target_definition.components.get(annos[i][5]).definition
            
            if i in indices:
                target_definition.sequenceAnnotations.remove(annos[i][2])

                if annos[i][5]:
                    target_definition.components.remove(annos[i][5])

                self.logger.debug('Removed %s at [%s, %s] in %s', feature_identity, annos[i][0], annos[i][1],
                              target_definition.identity)

                del annos[i]

    @classmethod
    def get_common_role(cls, roles):
        for role in roles:
            if role in cls.COMMON_ROLE_DICT:
                return cls.COMMON_ROLE_DICT[role]

        return 'sequence_feature'

    @classmethod
    def __ask_swap_annotations(cls, anno, sub_part_anno, target_definition):
        if anno[4] is None:
            anno_ID = anno[3]
        else:
            anno_ID = anno[4]

        if sub_part_anno[4] is None:
            sub_part_anno_ID = sub_part_anno[3]
        else:
            sub_part_anno_ID = sub_part_anno[4]

        select_message = 'Annotation {an} ({ai}) and\nsub-part annotation {sp} ({si})\nat [{st}, {en}] in {td} appear \
to be nearly identical.\nRemove second annotation and link first to sub-part? \
(0=no,1=yes):\n'.format(an=anno[2],
                        ai=anno_ID,
                        sp=sub_part_anno[2],
                        si=sub_part_anno_ID,
                        st=anno[0],
                        en=anno[1],
                        td=target_definition.identity)

        selected_message = input(select_message)

        try:
            return int(selected_message.strip()) == 1
        except ValueError:
            return False

    @classmethod
    def __select_annotations(cls, doc, target_definition, annos, ask_user=True, canonical_library=None, keep_flat=True):
        kept_indices = []

        feature_messages = []

        for i in range(0, len(annos)):
            if annos[i][5] is None:
                if annos[i][4] is None:
                    feature_ID = annos[i][3]
                else:
                    feature_ID = annos[i][4]

                feature_role = cls.get_common_role(annos[i][6])

                if ask_user:
                    if annos[i][7] and len(annos[i][7]) > 0:
                        if len(feature_role) > 0:
                            feature_messages.append('{nx}: {id} ({fi}, {ro}) at [{st}, {en}]. {de}'.format(nx=str(i), id=annos[i][2],
                                fi=feature_ID, ro=feature_role, st=annos[i][0], en=annos[i][1], de=annos[i][7]))
                        else:
                            feature_messages.append('{nx}: {id} ({fi}) at [{st}, {en}]. {de}'.format(nx=str(i), id=annos[i][2],
                                fi=feature_ID, st=annos[i][0], en=annos[i][1], de=annos[i][7]))
                    elif len(feature_role) > 0:
                        feature_messages.append('{nx}: {id} ({fi}, {ro}) at [{st}, {en}]'.format(nx=str(i), id=annos[i][2],
                            fi=feature_ID, ro=feature_role, st=annos[i][0], en=annos[i][1]))
                    else:
                        feature_messages.append('{nx}: {id} ({fi}) at [{st}, {en}]'.format(nx=str(i), id=annos[i][2],
                            fi=feature_ID, st=annos[i][0], en=annos[i][1]))
                elif keep_flat:
                    kept_indices.append(i)
            else:
                feature_identity = target_definition.components.get(annos[i][5]).definition

                if ask_user:
                    feature_definition = doc.getComponentDefinition(feature_identity)
                    
                    if feature_definition.name is None:
                        feature_ID = feature_definition.displayId
                    else:
                        feature_ID = feature_definition.name

                    feature_role = cls.get_common_role(feature_definition.roles)

                    feature_description = feature_definition.description

                    if feature_description and len(feature_description) > 0:
                        if len(feature_role) > 0:
                            feature_messages.append('{nx}: {id} ({fi}, {ro}) at [{st}, {en}]. {de}'.format(nx=str(i), id=feature_identity,
                                fi=feature_ID, ro=feature_role, st=annos[i][0], en=annos[i][1], de=feature_description))
                        else:
                            feature_messages.append('{nx}: {id} ({fi}) at [{st}, {en}]. {de}'.format(nx=str(i), id=feature_identity,
                                fi=feature_ID, st=annos[i][0], en=annos[i][1], de=feature_description))
                    elif len(feature_role) > 0:
                        feature_messages.append('{nx}: {id} ({fi}, {ro}) at [{st}, {en}]'.format(nx=str(i), id=feature_identity,
                            fi=feature_ID, ro=feature_role, st=annos[i][0], en=annos[i][1]))
                    else:
                        feature_messages.append('{nx}: {id} ({fi}) at [{st}, {en}]'.format(nx=str(i), id=feature_identity,
                            fi=feature_ID, st=annos[i][0], en=annos[i][1]))
                elif not canonical_library or canonical_library.has_feature(feature_identity):
                    kept_indices.append(i)

        if ask_user:
            if target_definition.name is None:
                target_ID = target_definition.displayId
            else:
                target_ID = target_definition.name

            select_message = 'There appear to be redundant features in {pi}:\n{fm}\nPlease select which ones to \
remove if any (comma-separated list of indices, for example 0,2,5):\n'.format(fm='\n'.join(feature_messages),
                                                                              pi=target_ID)

            selected_message = input(select_message)

            try:
                selected_indices = [int(si.strip()) for si in selected_message.split(',')]
            except ValueError:
                selected_indices = []

            return set(selected_indices)
        else:
            return set(range(0, len(annos))).difference(set(kept_indices))

    def __filter_annotations(self, annos, target_definition):
        for i in range(len(annos) - 1, -1, -1):
            if annos[i][5] is None:
                feature_identity = annos[i][2]

                feature_roles = annos[i][6]
            else:
                feature_identity = target_definition.components.get(annos[i][5]).definition

                feature_roles = set(self.feature_library.get_definition(feature_identity).roles)

            if len(feature_roles.intersection(self.roles)) == 0:
                target_definition.sequenceAnnotations.remove(annos[i][2])

                if annos[i][5]:
                    target_definition.components.remove(annos[i][5])

                self.logger.debug('Removed %s at [%s, %s] in %s', feature_identity, annos[i][0], annos[i][1],
                              target_definition.identity)

                del annos[i]

    def __are_swappable(self, anno, sub_part_anno, target_definition):
        if not anno[5] and sub_part_anno[5] and anno[0] == sub_part_anno[0] and anno[1] == sub_part_anno[1]:
            feature_identity = target_definition.components.get(sub_part_anno[5]).definition

            feature_definition = self.feature_library.get_definition(feature_identity)

            return (not Feature.has_non_generic_role(anno[6]) or
                    anno[6] == set(feature_definition.roles))
        else:
            return False

    def __swap_annotations(self, anno, sub_part_anno, target_definition):
        seq_anno = target_definition.sequenceAnnotations.get(anno[2])

        seq_anno.roles = []
        seq_anno.component = sub_part_anno[5]

        target_definition.sequenceAnnotations.remove(sub_part_anno[2])

        self.logger.debug('Removed %s at [%s, %s] in %s', sub_part_anno[2], sub_part_anno[0], sub_part_anno[1],
                      target_definition.identity)
        self.logger.debug('Modified %s at [%s, %s] in %s to refer to %s', anno[2], anno[0], anno[1],
                      target_definition.identity, sub_part_anno[5])

    # @classmethod
    # def __get_annotations(cls, doc, comp_definition):
    #     cut_annos = [(sa.locations.getCut().at, sa.locations.getCut().at, sa.identity, sa.displayId, sa.name, sa.component, set(sa.roles)) for sa in comp_definition.sequenceAnnotations if len(sa.locations) == 1 and sa.locations[0].getTypeURI() == SBOL_CUT]
    #     annos = [(sa.locations.getRange().start, sa.locations.getRange().end, sa.identity, sa.displayId, sa.name, sa.component, set(sa.roles)) for sa in comp_definition.sequenceAnnotations if len(sa.locations) == 1 and sa.locations[0].getTypeURI() == SBOL_RANGE]
            
    #     annos.extend(cut_annos)

    #     for sub_comp in comp_definition.components:
    #         try:
    #             sub_definition = doc.getComponentDefinition(comp_definition.identity)
    #         except RuntimeError:
    #             sub_definition = None

    #         if sub_definition is not None:
    #             annos.extend(cls.__get_annotations(doc, sub_definition))

    #     return annos

    @classmethod
    def __get_flat_annotation_indices(cls, anno_group):
        flat_indices = []

        for i in range(0, len(anno_group)):
            if anno_group[i][5] is None:
                flat_indices.append(i)

        return flat_indices

    def clean(self, feature_library, annotated_features, annotating_features):
        self.logger.info('Cleaning up')

        sub_definitions = set()

        for annotated_feature in annotated_features:
            annotated_doc = feature_library.get_document(annotated_feature.identity)
            annotated_definition = annotated_doc.getComponentDefinition(annotated_feature.identity)

            sub_IDs = set()
            temp_sub_definitions = set()

            for comp in annotated_definition.components:
                sub_IDs.add(comp.displayId)
                temp_sub_definitions.add(comp.definition)
            for seq_anno in annotated_definition.sequenceAnnotations:
                sub_IDs.add(seq_anno.displayId)

            parent_sub_IDs = set()

            for parent_identity in annotated_definition.wasDerivedFrom:
                parent_doc = feature_library.get_document(parent_identity)

                if parent_doc:
                    parent_definition = parent_doc.getComponentDefinition(parent_identity)

                    for comp in parent_definition.components:
                        parent_sub_IDs.add(comp.displayId)
                    for seq_anno in parent_definition.sequenceAnnotations:
                        parent_sub_IDs.add(seq_anno.displayId)

            if len(sub_IDs) == len(parent_sub_IDs) and len(sub_IDs) == len(sub_IDs.intersection(parent_sub_IDs)):
                annotated_doc.componentDefinitions.remove(annotated_feature.identity)

                self.logger.debug('Removed %s from %s', annotated_feature.identity, annotated_doc.name)
            else:
                sub_definitions.update(temp_sub_definitions)

        for annotating_feature in annotating_features:
            if annotating_feature.identity not in sub_definitions:
                annotating_doc = feature_library.get_document(annotating_feature.identity)
                annotating_definition = annotating_doc.getComponentDefinition(annotating_feature.identity)

                annotating_doc.componentDefinitions.remove(annotating_feature.identity)
                for seq_identity in annotating_definition.sequences:
                    try:
                        annotating_doc.sequences.remove(seq_identity)
                    except ValueError:
                        pass

                self.logger.debug('Removed %s from %s', annotating_feature.identity, annotating_doc.name)

        self.logger.info('Finished cleaning up')

    def prune(self, target_library, cover_offset, min_target_length, ask_user=True, canonical_library=None,
              delete_flat=False, keep_flat=True, target_features=[], auto_swap=False, require_sequence=True):
        target_identities = set()
        for target_feature in target_features:
            target_identities.add(target_feature.identity)

        for target in target_library.features:
            if ((not require_sequence or self.__has_min_length(target, min_target_length))
                    and (len(target_identities) == 0 or target.identity in target_identities)):
                self.logger.info('Pruning %s', target.identity)

                target_doc = target_library.get_document(target.identity)

                target_definition = target_doc.getComponentDefinition(target.identity)

                cut_annos = [(sa.locations.getCut().at, sa.locations.getCut().at, sa.identity, sa.displayId, sa.name, sa.component, set(sa.roles), sa.description) for sa in target_definition.sequenceAnnotations if len(sa.locations) == 1 and sa.locations[0].getTypeURI() == sbol2.SBOL_CUT]
                annos = [(sa.locations.getRange().start, sa.locations.getRange().end, sa.identity, sa.displayId, sa.name, sa.component, set(sa.roles), sa.description) for sa in target_definition.sequenceAnnotations if len(sa.locations) == 1 and sa.locations[0].getTypeURI() == sbol2.SBOL_RANGE]
                
                annos.extend(cut_annos)

                if len(self.roles) > 0:
                    self.__filter_annotations(annos, target_definition)

                annos.sort()

                grouped_annos = [[]]

                for anno in annos:
                    if len(grouped_annos) > 1 and self.__is_covered(anno, grouped_annos[-2], cover_offset):
                        grouped_annos[-2].append(anno)
                    elif self.__is_covered(anno, grouped_annos[-1], cover_offset):
                        grouped_annos[-1].append(anno)
                    else:
                        grouped_annos.append([anno])

                if delete_flat:
                    for anno_group in grouped_annos:
                        flat_indices = self.__get_flat_annotation_indices(anno_group)

                        self.__remove_annotations(flat_indices, anno_group, target_definition)

                for anno_group in grouped_annos:
                    if len(anno_group) > 1:
                        selected_indices = self.__select_annotations(target_doc, target_definition, anno_group,
                            ask_user, canonical_library, keep_flat)

                        self.__remove_annotations(selected_indices, anno_group, target_definition)

                for anno_group in grouped_annos:
                    if ask_user or auto_swap:
                        if len(anno_group) == 2:
                            if (self.__are_swappable(anno_group[0], anno_group[1], target_definition) and
                                    (auto_swap or
                                    self.__ask_swap_annotations(anno_group[0], anno_group[1], target_definition))):
                                self.__swap_annotations(anno_group[0], anno_group[1], target_definition)
                            elif (self.__are_swappable(anno_group[1], anno_group[0], target_definition) and
                                    (auto_swap or
                                    self.__ask_swap_annotations(anno_group[1], anno_group[0], target_definition))):
                                self.__swap_annotations(anno_group[1], anno_group[0], target_definition)

                for anno_group in grouped_annos:
                    if len(anno_group) > 1:
                        redundant_URIs = []

                        for anno in anno_group:
                            if anno[5]:
                                redundant_URIs.append(target_definition.components.get(anno[5]).definition)
                            else:
                                redundant_URIs.append(anno[2])

                        self.logger.debug('Detected potentially redundant sub-parts in %s: %s.',
                                          target_definition.identity, 
                                          str(redundant_URIs))

                self.logger.info('Finished pruning %s', target.identity)
    