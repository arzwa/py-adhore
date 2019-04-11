import logging
import os


def default_adhore_conf(**kwargs):
    default_conf = {
        "gap_size": 30,
        "q_value": 0.75,
        "cluster_gap": 35,
        "prob_cutoff": 0.01,
        "anchor_points": 3,
        "alignment_method": "gg2",
        "level_2_only": "false",
        "table_type": "family",
        "multiple_hypothesis_correction": "FDR",
        "visualizeGHM": "false",
        "visualizeAlignment": "true"
    }
    for k, v in kwargs.items():
        if k in default_conf:
            default_conf[k] = v
        else:
            logging.warning("{} is not a vali I-ADHoRe parameter, ignoring")
    return default_conf


def write_adhore_config(conf, path):
    """
    Write the I-ADHoRe config file.
    """
    with open(path, "w") as f:
        for k, v in conf.items():
            if k == "lists":
                for genome in v:
                    f.write("genome={}\n".format(genome["genome"]))
                    f.write("\n".join(genome["lists"]))
                    f.write("\n\n")
            else:
                f.write("{}={}\n".format(k, v))
    return os.path.abspath(path)


def parse_feat_attr(features, attributes, species):
    """
    Parse and interpret the features and attributes option settings for CLI
    programs.
    """
    if len(features.split(",")) == 1:
        logging.warning("# features different from # species, will use"
        " {} for all species".format(features))
        feat = [features.split(";") for i in species]
    else:
        feat = [x.split(";") for x in features.split(",")]
    if len(attributes.split(",")) == 1:
        logging.warning("# attributes different from # species, will use"
        " {} for all species".format(attributes))
        attr = [attributes for i in species]
    else:
        attr = attributes.split(",")
    logging.info("Species and their gff specs:")
    for i, s in enumerate(species):
        logging.info(" .. {}: feat. = {}; attr. = {}".format(
            s, feat[i], attr[i]))
    return feat, attr
