#!/usr/bin/env python3
"""Render template.xml from YAML (constants, groups, fork cases; anchors supported)."""
import argparse
import yaml
import itertools
import copy
from pathlib import Path
from typing import Any, Dict, List, Tuple

_C_RESET = "\033[0m"
_C_GREEN = "\033[32m"
_C_RED = "\033[31m"

def _green(s: str) -> str:
    return f"{_C_GREEN}{s}{_C_RESET}"

def _red(s: str) -> str:
    return f"{_C_RED}{s}{_C_RESET}"


DEFAULT_CONFIG: Dict[str, Any] = {
    "Simulation": {
        "SimType": "Standard",
        "Precision": "double",
    },
    "DomainOrigin": {"X": 0.0, "Y": 0.0, "Z": 0.0},
    "DomainSize": {"X": 1.0, "Y": 1.0, "Z": 1.0},

    "CollisionDetection": {
        "NeighborListType": "LinkedCell",
        "LinkedCellType": "Host",
        "CellSizeFactor": 1.0,
        "UpdateFrequency": 0,
        "SortingFrequency": 0,
        "BoundingVolumeType": "OBB",
        "NarrowPhaseType": "GJK",
        "EnableTimings": "false",
        "UseRelativeTransformations": "false",
        "NarrowPhaseAcceleration": "false"
    },

    "ContactModel": {
        "ContactModelType": "Hooke",
        "MaterialA": "MatP",
        "MaterialB": "MatP",
        "Parameters": 'kn="1.2e8" en="0.8" etat="100.0" muc="0.5" kr="0.0"'
    },
    "ContactModelsList": [],

    "ContactForceModels": {
        "EnableTimings": "false",
        "Compaction": "false"
    },

    "Temporal": {"TStart": 0.0, "TEnd": 1.0, "DT": 1e-3, "TimeIntegrationType": "FirstOrderExplicit"},
    
    "Particles": {
        "NumParticles": 1,
        "ParticleDensity": 1000.0,
        "ParticleMaterial": "MatP",
        "CrustThickness": 0.0002,
        "ParticleShape": '<Sphere Radius="0.05"/>' ,
        "AngularPosition": {"aX": 0.0, "aY": 0.0, "aZ": 0.0}
    },
    "ParticlesList": [],
    
    "Obstacle": {
        "ObstacleMaterial": "rigid",
        "ObstacleShape": '<Vertex X="0" Y="0" Z="0"/>' ,
        "ObstaclePosition": {"X": 0.0, "Y": 0.0, "Z": 0.0},
        "ObstacleAngularPosition": {"aX": 0.0, "aY": 0.0, "aZ": 0.0}
    },
    "ObstaclesList": [],

    "Forces": {"Gravity": {"GX": 0.0, "GY": 0.0, "GZ": -9.81}},

    "SimulationSettings": {
        "ForceInsertion": "0",
        "VerbosityFrequency": 1000,
        "TimingsEnable": "false",
        "OutputPrecision": 6,
        "PositionSection": "",
        "OrientationSection": "",
        "VelocitySection": "",
        "AngularVelocitySection": "",
        "WritersSection": ""
    },

    "PostProcessing": {
        "TimeSave": {"Start": 0.0, "End": 1.0, "DT": 1e-3}
    },
    "Insertion": {
        "ForceInsertion": 0,
        "InitialPosition": {"Type": "Random", "Seed": {"Mode": "UserDefined", "Value": 1}},
        "InitialOrientation": {"Type": "Random", "Seed": {"Mode": "UserDefined", "Value": 1}},
        "InitialVelocity": {"Type": "Zero"},
        "InitialAngularVelocity": {"Type": "Zero"}
    }
}


def deep_merge(a: Dict[str, Any], b: Dict[str, Any]) -> Dict[str, Any]:
    """Recursive merge of b into a copy of a (b wins on conflicts)."""
    out = copy.deepcopy(a)
    for k, v in b.items():
        if k in out and isinstance(out[k], dict) and isinstance(v, dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


class TemplatePopulator:
    def __init__(self, config: Dict[str, Any]):
        self.config = copy.deepcopy(config)
        # Merge scalar `Particles` / `Obstacle` into each list entry when both exist.
        if isinstance(self.config.get("Particles"), dict) and isinstance(self.config.get("ParticlesList"), list):
            part_base = self.config.pop("Particles")
            for i, p in enumerate(self.config.get("ParticlesList", [])):
                self.config["ParticlesList"][i] = deep_merge(part_base, p)
        if isinstance(self.config.get("Obstacle"), dict) and isinstance(self.config.get("ObstaclesList"), list):
            obs_base = self.config.pop("Obstacle")
            for i, o in enumerate(self.config.get("ObstaclesList", [])):
                self.config["ObstaclesList"][i] = deep_merge(obs_base, o)

    def _flatten_for_template(self) -> Dict[str, Any]:
        """Map nested config to template placeholder keys."""
        m: Dict[str, Any] = {}

        m["sim_type"] = self.config.get("Simulation", {}).get("SimType", "")
        m["precision"] = self.config.get("Simulation", {}).get("Precision", "")

        for axis in ("X", "Y", "Z"):
            m[f"domain_origin_{axis}"] = self.config.get("DomainOrigin", {}).get(axis, "")
            m[f"domain_size_{axis}"] = self.config.get("DomainSize", {}).get(axis, "")

        cd = self.config.get("CollisionDetection", {})
        m["neighbor_list_type"] = cd.get("NeighborListType", "")
        m["linked_cell_type"] = cd.get("LinkedCellType", "")
        m["cell_size_factor"] = cd.get("CellSizeFactor", "")
        m["update_frequency"] = cd.get("UpdateFrequency", "")
        m["sorting_frequency"] = cd.get("SortingFrequency", "")
        m["bounding_volume_type"] = cd.get("BoundingVolumeType", "")
        m["narrow_phase_type"] = cd.get("NarrowPhaseType", "")
        m["cd_enable_timings"] = cd.get("EnableTimings", "false")
        m["cd_use_relative_transformations"] = cd.get("UseRelativeTransformations", "false")
        m["narrow_phase_acceleration"] = cd.get("NarrowPhaseAcceleration", "false")

        cfm_cfg = self.config.get("ContactForceModels", {})
        m["cfm_enable_timings"] = cfm_cfg.get("EnableTimings", "false")
        m["cfm_compaction"] = cfm_cfg.get("Compaction", "false")

        cm = self.config.get("ContactModel", {})
        m["contact_model_type"] = cm.get("ContactModelType", "")
        m["contact_model_materialA"] = cm.get("MaterialA", "")
        m["contact_model_materialB"] = cm.get("MaterialB", "")
        m["contact_model_parameters"] = cm.get("Parameters", "")

        t = self.config.get("Temporal", {})
        m["t_start"] = t.get("TStart", "")
        m["t_end"] = t.get("TEnd", "")
        m["dt"] = t.get("DT", "")
        m["time_integration_type"] = t.get("TimeIntegrationType", "")

        p = self.config.get("Particles", {})
        m["num_particles"] = p.get("NumParticles", "")
        m["particle_density"] = p.get("ParticleDensity", "")
        m["particle_material"] = p.get("ParticleMaterial", "")
        m["crust_thickness"] = p.get("CrustThickness", "")
        m["particle_shape"] = p.get("ParticleShape", "")
        for ax in ("aX", "aY", "aZ"):
            m[f"particle_angular_position_{ax}"] = p.get("AngularPosition", {}).get(ax, "")

        o = self.config.get("Obstacle", {})
        m["obstacle_material"] = o.get("ObstacleMaterial", "")
        m["obstacle_shape"] = o.get("ObstacleShape", "")
        for axis in ("X", "Y", "Z"):
            m[f"obstacle_position_{axis}"] = o.get("ObstaclePosition", {}).get(axis, "")
        for ax in ("aX", "aY", "aZ"):
            m[f"obstacle_angular_position_{ax}"] = o.get("ObstacleAngularPosition", {}).get(ax, "")

        g = self.config.get("Forces", {}).get("Gravity", {})
        m["gravity_GX"] = g.get("GX", "")
        m["gravity_GY"] = g.get("GY", "")
        m["gravity_GZ"] = g.get("GZ", "")

        ss = self.config.get("SimulationSettings", {})
        m["force_insertion"] = ss.get("ForceInsertion", "false")
        m["verbosity_frequency"] = ss.get("VerbosityFrequency", 1000)
        m["simulation_timings_enable"] = ss.get("TimingsEnable", "false")

        def build_section(key: str, tag: str) -> str:
            val = ss.get(key, "")
            if isinstance(val, list):
                if val and isinstance(val[0], dict):
                    lines = self._build_writers_lines(val)
                    return "\n".join(lines)
                return "\n" + "\n".join(val)
            if not val and tag == "Writer":
                return "\n" + "            "
            return str(val or "")

        m["position_section"] = build_section("PositionSection", "Position")
        m["orientation_section"] = build_section("OrientationSection", "Orientation")
        m["velocity_section"] = build_section("VelocitySection", "Velocity")
        m["angular_velocity_section"] = build_section("AngularVelocitySection", "AngularVelocity")
        m["writers_section"] = build_section("WritersSection", "Writer")

        pp_root = self.config.get("PostProcessing", {})
        pp = pp_root.get("TimeSave", {})
        m["time_save_start"] = pp.get("Start", "")
        m["time_save_end"] = pp.get("End", "")
        m["time_save_dt"] = pp.get("DT", "")

        return m

    def _build_writers_lines(self, writers_list) -> list:
        """Writer XML lines (RawData / Paraview) from structured dict list."""
        lines = []
        default_precision = self.config.get("SimulationSettings", {}).get("OutputPrecision", 6)
        for w in writers_list:
            typ = w.get("Type", "Raw")
            directory = w.get("Directory", "./results")
            root = w.get("RootName", "output")
            if typ.lower() in ("raw", "rawdata"):
                precision = w.get("Precision", default_precision)
                lines.append(f'            <RawData Directory ="{directory}" RootName="{root}" Precision="{precision}"/>')
            elif typ.lower() in ("paraview", "para", "pv"):
                lines.append(f'            <Paraview Directory ="{directory}" RootName="{root}"/>')
            else:
                lines.append(f'            <RawData Directory ="{directory}" RootName="{root}" Precision="{default_precision}"/>')
        return lines

    def render_from_template(self, tpl: str) -> str:
        mapping = self._flatten_for_template()

        class SafeDict(dict):
            def __missing__(self, key):
                return ""

        try:
            rendered = tpl.format_map(SafeDict(mapping))
        except Exception as e:
            raise RuntimeError(f"Failed to render template: {e}")

        if not self.config.get("ParticlesList"):
            rendered = self._remove_block(rendered, "Particles")
        if not self.config.get("ObstaclesList"):
            rendered = self._remove_block(rendered, "Obstacles")

        rendered = self._replace_block_if_list(rendered, "Particles", self.config.get("ParticlesList"), self._build_particles_block)
        rendered = self._replace_block_if_list(rendered, "Obstacles", self.config.get("ObstaclesList"), self._build_obstacles_block)
        rendered = self._replace_block_if_list(rendered, "ContactForceModels", self.config.get("ContactModelsList"), self._build_contact_models_block)

        rendered = self._replace_insertion_block(rendered, self.config.get("Insertion"))

        return rendered

    def _replace_block_if_list(self, rendered: str, tag: str, lst, builder_func) -> str:
        if not lst:
            return rendered
        block = builder_func(lst)
        # Prefix match: opening tag may include attributes.
        start_tag_prefix = f"<{tag}"
        end_tag = f"</{tag}>"
        start_idx = rendered.find(start_tag_prefix)
        if start_idx != -1:
            open_end_idx = rendered.find(">", start_idx)
            end_idx = rendered.find(end_tag, open_end_idx) if open_end_idx != -1 else -1
        else:
            end_idx = -1
        if start_idx != -1 and end_idx != -1:
            line_start = rendered.rfind("\n", 0, start_idx) + 1
            indent = rendered[line_start:start_idx]
            indented_block = indent + block.replace("\n", "\n" + indent)
            return rendered[: line_start] + indented_block + rendered[end_idx + len(end_tag) :]
        return rendered + "\n" + block

    def _build_particles_block(self, particles_list) -> str:
        """Full <Particles>...</Particles> from a list of particle dicts."""
        parts = ["<Particles>"]
        for p in particles_list:
            attrs = []
            if "Number" in p:
                attrs.append(f'Number="{p["Number"]}"')
            if "Density" in p:
                attrs.append(f'Density="{p["Density"]}"')
            if "Material" in p:
                attrs.append(f'Material="{p["Material"]}"')
            attr_str = " " + " ".join(attrs) if attrs else ""

            # Use PascalCase fields only
            crust = p.get("CrustThickness", 0.0)
            shape = p.get("ParticleShape", "")
            ang = p.get("AngularPosition", {})
            aX = ang.get("aX", "0.0")
            aY = ang.get("aY", "0.0")
            aZ = ang.get("aZ", "0.0")

            part = []
            part.append(f"    <Particle{attr_str}>")
            part.append(f"        <Convex CrustThickness=\"{crust}\">")
            if shape:
                for line in str(shape).splitlines():
                    part.append("            " + line)
            part.append("        </Convex>")

            pos = p.get("Position") or p.get("ParticlePosition") or {}
            X = pos.get("X", 0.0)
            Y = pos.get("Y", 0.0)
            Z = pos.get("Z", 0.0)
            part.append("        <Transformation>")
            part.append(f"            <Centre X=\"{X}\" Y=\"{Y}\" Z=\"{Z}\"/>")
            part.append(f"            <AngularPosition Type=\"Angles\" aX=\"{aX}\" aY=\"{aY}\" aZ=\"{aZ}\"/>")
            part.append("        </Transformation>")
            part.append("    </Particle>")
            parts.extend(part)

        parts.append("</Particles>")
        return "\n".join(parts)

    def _build_obstacles_block(self, obstacles_list) -> str:
        """Full <Obstacles>...</Obstacles> from obstacle dicts."""
        parts = ["<Obstacles>"]
        for o in obstacles_list:
            mat = o.get("Material") or o.get("ObstacleMaterial", "")
            crust = o.get("CrustThickness", 0.0)
            shape = o.get("ObstacleShape", "")
            pos = o.get("ObstaclePosition") or o.get("Position", {})
            ang = o.get("ObstacleAngularPosition") or o.get("AngularPosition", {})
            aX = ang.get("aX", "0.0")
            aY = ang.get("aY", "0.0")
            aZ = ang.get("aZ", "0.0")
            X = pos.get("X", "0.0")
            Y = pos.get("Y", "0.0")
            Z = pos.get("Z", "0.0")

            parts.append(f"    <Obstacle Material=\"{mat}\">")
            parts.append(f"        <Convex CrustThickness=\"{crust}\">")
            if shape:
                for line in str(shape).splitlines():
                    parts.append("            " + line)
            parts.append("        </Convex>")
            parts.append("        <Transformation>")
            parts.append(f"            <Centre X=\"{X}\" Y=\"{Y}\" Z=\"{Z}\"/>")
            parts.append(f"            <AngularPosition Type=\"Angles\" aX=\"{aX}\" aY=\"{aY}\" aZ=\"{aZ}\"/>")
            parts.append("        </Transformation>")
            parts.append("    </Obstacle>")

        parts.append("</Obstacles>")
        return "\n".join(parts)

    def _build_contact_models_block(self, models_list) -> str:
        """Full <ContactForceModels>...</ContactForceModels> from model dicts."""
        cfm_cfg = self.config.get("ContactForceModels", {})
        enable_timings = cfm_cfg.get("EnableTimings", "false")
        compaction = cfm_cfg.get("Compaction", "false")
        parts = [f'<ContactForceModels EnableTimings="{enable_timings}" Compaction="{compaction}">']
        for m in models_list:
            typ = m.get("Type", m.get("ContactModelType", ""))
            matA = m.get("MaterialA", "")
            matB = m.get("MaterialB", "")
            params = m.get("Parameters", "")
            parts.append(f"    <ContactForceModel Type=\"{typ}\">")
            parts.append(f"        <Material materialA=\"{matA}\" materialB=\"{matB}\"/>")
            if params:
                parts.append(f"        <Parameters {params}/>")
            parts.append("    </ContactForceModel>")
        parts.append("</ContactForceModels>")
        return "\n".join(parts)

    def _replace_insertion_block(self, rendered: str, insertion_cfg) -> str:
        if not insertion_cfg:
            return rendered
        block = self._build_insertion_block(insertion_cfg)
        start_tag = "<ParticleInsertion"
        end_tag = "</ParticleInsertion>"
        start_idx = rendered.find(start_tag)
        end_idx = rendered.find(end_tag, start_idx)
        if start_idx != -1 and end_idx != -1:
            line_start = rendered.rfind("\n", 0, start_idx) + 1
            indent = rendered[line_start:start_idx]
            indented_block = indent + block.replace("\n", "\n" + indent)
            return rendered[: line_start] + indented_block + rendered[end_idx + len(end_tag) :]
        return rendered + "\n" + block

    def _remove_block(self, rendered: str, tag: str) -> str:
        """Strip `<tag>...</tag>` (whole lines); no-op if not found."""
        start_tag = f"<{tag}>"
        end_tag = f"</{tag}>"
        start_idx = rendered.find(start_tag)
        end_idx = rendered.find(end_tag, start_idx)
        if start_idx != -1 and end_idx != -1:
            line_start = rendered.rfind("\n", 0, start_idx) + 1
            return rendered[:line_start] + rendered[end_idx + len(end_tag) :]
        return rendered

    def _build_insertion_block(self, cfg) -> str:
        """<ParticleInsertion> from cfg (Random / File / Constant / Zero)."""
        fi = cfg.get("ForceInsertion", False)
        fi_attr = "1" if fi else "0"
        parts = [f"<ParticleInsertion ForceInsertion=\"{fi_attr}\">"]

        def build_seed_attr(seed):
            if not seed:
                return "", {}
            mode = seed.get("Mode", "Default")
            attrs = {"Seed": mode}
            if mode == "UserDefined":
                attrs["Value"] = seed.get("Value", 1)
            return ' '.join(f'{k}="{v}"' for k, v in attrs.items()), attrs

        def build_windows(windows):
            if not windows:
                domain = self.config.get("DomainSize", {})
                maxx = domain.get("X", 1.0)
                maxy = domain.get("Y", maxx)
                maxz = domain.get("Z", maxx)
                return [
                    {
                        "Type": "Box",
                        "MinPoint": {"X": 0.0, "Y": 0.0, "Z": 0.0},
                        "MaxPoint": {"X": maxx, "Y": maxy, "Z": maxz},
                    }
                ]
            return windows

        for key in ("InitialPosition", "InitialOrientation", "InitialVelocity", "InitialAngularVelocity"):
            node = cfg.get(key, {"Type": "Zero"})
            ntype = node.get("Type", "Zero")
            if ntype == "Random":
                seed_str, _ = build_seed_attr(node.get("Seed"))
                attr = f'Type="Random"'
                if seed_str:
                    attr = f'{attr} {seed_str}'
                parts.append(f"    <{key} {attr}>")
                parts.append("        <Windows>")
                for w in build_windows(node.get("Windows")):
                    wtype = w.get("Type", "Box")
                    parts.append(f"            <Window Type=\"{wtype}\">")
                    if wtype == "Box":
                        mp = w.get("MinPoint", {})
                        xp = w.get("MaxPoint", {})
                        parts.append(f"                <MinPoint X=\"{mp.get('X',0)}\" Y=\"{mp.get('Y',0)}\" Z=\"{mp.get('Z',0)}\"/>")
                        parts.append(f"                <MaxPoint X=\"{xp.get('X',0)}\" Y=\"{xp.get('Y',0)}\" Z=\"{xp.get('Z',0)}\"/>")
                    else:
                        for k, v in w.items():
                            if k == "Type":
                                continue
                            parts.append(f"                <{k}>{v}</{k}>")
                    parts.append(f"            </Window>")
                parts.append("        </Windows>")
                parts.append(f"    </{key}>")
            elif ntype == "File":
                name = node.get("Name", "")
                parts.append(f"    <{key} Type=\"File\" Name=\"{name}\"/>")
            elif ntype == "Constant":
                x = node.get("X") if node.get("X") is not None else node.get("Value", {}).get("X", 0)
                y = node.get("Y") if node.get("Y") is not None else node.get("Value", {}).get("Y", 0)
                z = node.get("Z") if node.get("Z") is not None else node.get("Value", {}).get("Z", 0)
                parts.append(f"    <{key} Type=\"Constant\" X=\"{x}\" Y=\"{y}\" Z=\"{z}\"/>")
            else:
                parts.append(f"    <{key} Type=\"Zero\"/>")

        parts.append("</ParticleInsertion>")
        return "\n".join(parts)

def build_configs_from_yaml(data: Dict[str, Any]) -> List[Tuple[str, Dict[str, Any]]]:
    """(output_basename, merged config) pairs from YAML Constant/Groups/Fork."""
    const = data.get("Constant", {})
    groups = data.get("Groups", {}) or {}
    fork = data.get("Fork", [])

    results = []
    if not fork:
        cfg = deep_merge(DEFAULT_CONFIG, const)
        results.append(("default", cfg))
        return results

    for idx, entry in enumerate(fork):
        name_base = entry.get("name") or f"case_{idx}"
        overrides = entry.get("overrides", {}) or {}

        grp_spec = entry.get("groups") or entry.get("use_groups") or []

        def _normalize_cfg(cfg: Dict[str, Any]):
            cd = cfg.get("CollisionDetection", {})
            if "SortFrequency" in cd and "SortingFrequency" not in cd:
                cd["SortingFrequency"] = cd.pop("SortFrequency")
            cfg["CollisionDetection"] = cd
            return cfg

        def _find_list_leaves(d, path=None):
            path = path or []
            leaves = []
            if isinstance(d, dict):
                for k, v in d.items():
                    leaves.extend(_find_list_leaves(v, path + [k]))
            else:
                if isinstance(d, list) and d and not any(isinstance(x, dict) for x in d):
                    leaves.append((path, d))
            return leaves

        def _set_in_dict(d, path, val):
            if not path:
                return val
            key = path[0]
            if len(path) == 1:
                d[key] = val
                return d
            if key not in d or not isinstance(d[key], dict):
                d[key] = {}
            d[key] = _set_in_dict(d.get(key, {}), path[1:], val)
            return d

        def _expand_range_shorthand_inplace(obj):
            if isinstance(obj, dict):
                for k, v in list(obj.items()):
                    if isinstance(v, dict) and "range" in v and isinstance(v["range"], list):
                        obj[k] = v["range"]
                    else:
                        _expand_range_shorthand_inplace(v)
            elif isinstance(obj, list):
                for item in obj:
                    _expand_range_shorthand_inplace(item)

        def _apply_overrides(base_cfg, overrides_dict):
            c = copy.deepcopy(base_cfg)
            ov = copy.deepcopy(overrides_dict)

            def _merge_into_transformation(target_elem, elem_override):
                if not isinstance(elem_override, dict):
                    return target_elem
                ov_copy = copy.deepcopy(elem_override)
                if "Transformation" in ov_copy:
                    transf = ov_copy.pop("Transformation")
                    target_elem["Transformation"] = deep_merge(target_elem.get("Transformation", {}), transf)

                transf_keys = ["AngularPosition", "Centre", "CentrePosition"]
                transf_collect = {}
                for k in transf_keys:
                    if k in ov_copy:
                        transf_collect[k] = ov_copy.pop(k)
                if transf_collect:
                    target_elem["Transformation"] = deep_merge(target_elem.get("Transformation", {}), transf_collect)
                    if "AngularPosition" in transf_collect:
                        target_elem["AngularPosition"] = deep_merge(target_elem.get("AngularPosition", {}), transf_collect["AngularPosition"])
                    if "Centre" in transf_collect:
                        pos = transf_collect.get("Centre", {})
                        target_elem["Position"] = deep_merge(target_elem.get("Position", {}), pos)

                return deep_merge(target_elem, ov_copy)

            if "Particles" in ov and isinstance(c.get("ParticlesList"), list):
                part_ov = ov.pop("Particles")
                for i, p in enumerate(c["ParticlesList"]):
                    merged = _merge_into_transformation(copy.deepcopy(p), part_ov)
                    c["ParticlesList"][i] = merged

            if "Obstacles" in ov and isinstance(c.get("ObstaclesList"), list):
                obs_ov = ov.pop("Obstacles")
                for i, o in enumerate(c["ObstaclesList"]):
                    merged = _merge_into_transformation(copy.deepcopy(o), obs_ov)
                    c["ObstaclesList"][i] = merged

            return deep_merge(c, ov)

        if not grp_spec:
            cfg = deep_merge(DEFAULT_CONFIG, const)
            cfg = _apply_overrides(cfg, overrides)
            cfg = _normalize_cfg(cfg)
            results.append((name_base, cfg))
            continue

        if all(not isinstance(g, list) for g in grp_spec):
            for g in grp_spec:
                safe_g = str(g).replace(" ", "_")
                name = name_base + "_" + safe_g
                cfg = deep_merge(DEFAULT_CONFIG, const)
                cfg = deep_merge(cfg, groups.get(g, {}))

                ov_copy = copy.deepcopy(overrides)
                _expand_range_shorthand_inplace(ov_copy)
                list_leaves = _find_list_leaves(ov_copy)
                if not list_leaves:
                    cfg2 = _apply_overrides(cfg, ov_copy)
                    cfg2 = _normalize_cfg(cfg2)
                    results.append((name, cfg2))
                else:
                    paths, lists = zip(*list_leaves)
                    for values in itertools.product(*lists):
                        ov = copy.deepcopy(ov_copy)
                        suffix_parts = []
                        for p, v in zip(paths, values):
                            _set_in_dict(ov, p, v)
                            suffix_parts.append("".join([str(x) for x in p]) + "=" + str(v))
                        name2 = name + "_" + "+".join(suffix_parts)
                        cfg2 = _apply_overrides(cfg, ov)
                        cfg2 = _normalize_cfg(cfg2)
                        results.append((name2, cfg2))
            continue

        group_options = []
        for g in grp_spec:
            if isinstance(g, list):
                group_options.append(g)
            else:
                group_options.append([g])

        for combo in itertools.product(*group_options):
            safe_combo = [str(x).replace(" ", "_") for x in combo]
            name = name_base + "_" + "-".join(safe_combo)
            cfg = deep_merge(DEFAULT_CONFIG, const)
            for g in combo:
                cfg = deep_merge(cfg, groups.get(g, {}))

            ov_copy = copy.deepcopy(overrides)
            _expand_range_shorthand_inplace(ov_copy)
            list_leaves = _find_list_leaves(ov_copy)
            if not list_leaves:
                cfg2 = _apply_overrides(cfg, ov_copy)
                cfg2 = _normalize_cfg(cfg2)
                results.append((name, cfg2))
            else:
                paths, lists = zip(*list_leaves)
                for values in itertools.product(*lists):
                    ov = copy.deepcopy(ov_copy)
                    suffix_parts = []
                    for p, v in zip(paths, values):
                        _set_in_dict(ov, p, v)
                        suffix_parts.append("".join([str(x) for x in p]) + "=" + str(v))
                    name2 = name + "_" + "+".join(suffix_parts)
                    cfg2 = _apply_overrides(cfg, ov)
                    cfg2 = _normalize_cfg(cfg2)
                    results.append((name2, cfg2))
    return results


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate XML files from YAML configuration",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  python generate_tests.py -i config.yaml
  python generate_tests.py -i tests.yaml -t custom_template.xml -d output/
"""
    )
    parser.add_argument("--input", "-i", type=Path, help="YAML file with test configuration", required=False)
    parser.add_argument("--template", "-t", type=Path, default=Path(__file__).parent / "template.xml",
                        help="Template XML path (default: template.xml)")
    parser.add_argument("--out-dir", "-d", type=Path, default=Path(__file__).parent / "generated",
                        help="Output directory for generated XML files (default: generated/)")

    args = parser.parse_args()

    if not args.input:
        data = {}
    else:
        with open(args.input, 'r') as f:
            data = yaml.safe_load(f)
        if data is None:
            data = {}

        cases = build_configs_from_yaml(data)
    if args.input:
        stem = args.input.stem
        if len(cases) == 1:
            cases = [(stem, cases[0][1])]
        else:
            cases = [(f"{stem}_{name}", cfg) for (name, cfg) in cases]
    tpl = Path(args.template).read_text()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for name, case_cfg in cases:
        pop = TemplatePopulator(case_cfg)
        rendered = pop.render_from_template(tpl)
        out_path = out_dir / f"{name}.xml"
        out_path.write_text(rendered)
        print(f"Wrote: {_green(str(out_path))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
