#!/usr/bin/env python3
"""Render the WorksGood task graph to Mermaid diagram and table documentation.

Usage: python3 scripts/render_taskgraph.py
"""

import sys
import os
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))
from build_taskgraph import TASKS


def generate_mermaid() -> str:
    """Generate Mermaid flowchart syntax for the task dependencies."""
    lines = ["%% Task Graph for trajectory-replication research project", "graph TD"]

    # Group tasks by their first tag or "other"
    tag_groups = {}
    for task in TASKS:
        group = task.tags[0] if task.tags else "other"
        if group not in tag_groups:
            tag_groups[group] = []
        tag_groups[group].append(task)

    # Generate subgraph sections for each tag group
    for group, tasks in tag_groups.items():
        lines.append(f"    subgraph {group}")
        for task in tasks:
            # Create sanitized node name that can be used in Mermaid
            safe_id = task.id.replace("-", "_")
            lines.append(f"        {safe_id}[{task.title[:40]}]")
        lines.append(f"    end")

    # Add dependencies (edges)
    lines.append("")
    lines.append("    %% Dependencies")
    for task in TASKS:
        for dep in task.deps:
            safe_from = dep.replace("-", "_")
            safe_to = task.id.replace("-", "_")
            lines.append(f"    {safe_from} --> {safe_to}")

    return "\n".join(lines)


def generate_table() -> str:
    """Generate a Markdown table of all tasks with dependencies."""
    lines = [
        "| ID | Title | Dependencies | Tags | Hours | Deliverable |",
        "|----|-------|--------------|------|-------|-------------|"
    ]

    for task in TASKS:
        deps_str = ", ".join(f"`{dep}`" for dep in task.deps) if task.deps else "—"
        tags_str = ", ".join(f"`{tag}`" for tag in task.tags) if task.tags else "—"
        lines.append(f"| `{task.id}` | {task.title} | {deps_str} | {tags_str} | {task.hours} | `{task.deliverable}` |")

    return "\n".join(lines)


def generate_decisions() -> str:
    """Generate the DECISIONS.md file content."""
    decision_points = [
        {
            "number": 1,
            "question": "Can stable ancestry trajectories actually be inferred?",
            "outcome_if_no": "reformulate",
            "responsible_tasks": ["traj-validate"],
            "evidence_criteria": ["traj-validation.md shows stable structure", "traj-quality-report"]
        },
        {
            "number": 2,
            "question": "Under neutral demographic simulations, is the replication statistic calibrated?",
            "outcome_if_no": "fix before touching biological loci",
            "responsible_tasks": ["sim-sweeps", "sim-ctrl-null"],
            "evidence_criteria": ["neutral FPR ≤ 0.05", "calibration plots"]
        },
        {
            "number": 3,
            "question": "Does the method recover simulated repeated selection better than Akbari without greatly increasing false positives?",
            "outcome_if_no": "the method has not solved the problem",
            "responsible_tasks": ["sim-sweeps"],
            "evidence_criteria": ["power comparison", "ROC curves", "effect-size bias"]
        },
        {
            "number": 4,
            "question": "Do AMY and canonical loci behave sensibly?",
            "outcome_if_no": "investigate before scaling",
            "responsible_tasks": ["locus-amy", "locus-canonical"],
            "evidence_criteria": ["AMY trajectory patterns", "known locus replication"]
        },
        {
            "number": 5,
            "question": "Does the Akbari reanalysis reveal a meaningful low-Akbari/high-replication class?",
            "outcome_if_no": "conclude the hypothesised over-correction appears limited under this test",
            "responsible_tasks": ["locus-akbari"],
            "evidence_criteria": ["class exists and replicates", "enrichment statistics"]
        }
    ]

    lines = [
        "# Decision Points",
        "",
        "Critical go/no-go decisions for the trajectory-replication research project.",
        "",
        "| Decision | Question | If No | Responsible Tasks | Evidence Required | Status |",
        "|----------|----------|-------|-------------------|-------------------|--------|"
    ]

    for dp in decision_points:
        tasks_str = ", ".join(f"`{t}`" for t in dp["responsible_tasks"])
        evidence_str = ", ".join(f"{c}" for c in dp["evidence_criteria"])
        lines.append(f"| **{dp['number']}** | {dp['question']} | {dp['outcome_if_no']} | {tasks_str} | {evidence_str} | PENDING |")

    lines.extend([
        "",
        "## Decision Tracking",
        "",
        "Each decision will be recorded here with:",
        "- **Status**: PENDING / IN_PROGRESS / GO / NO_GO",
        "- **Date**: When the decision was made",
        "- **Evidence**: Link(s) to the supporting analysis/documentation",
        "- **Rationale**: Brief explanation of the reasoning",
        "",
        "This file should be updated as each decision point is reached and resolved."
    ])

    return "\n".join(lines)


def main() -> int:
    """Generate both TASKGRAPH.md and DECISIONS.md"""

    # Generate TASKGRAPH.md
    taskgraph_content = [
        "# Task Graph",
        "",
        "Dependency graph and task listing for the trajectory-replication research project.",
        "",
        "## Mermaid Flowchart",
        "",
        "```mermaid",
        generate_mermaid(),
        "```",
        "",
        "## Task Table",
        "",
        generate_table(),
        "",
        "## Legend",
        "",
        "- **ID**: Task identifier used in WG commands",
        "- **Title**: Human-readable description",
        "- **Dependencies**: Tasks that must complete before this task can start",
        "- **Tags**: Categorization labels",
        "- **Hours**: Estimated effort (person-hours)",
        "- **Deliverable**: Expected output file/path",
        "",
        "*Generated by `scripts/render_taskgraph.py`*"
    ]

    Path("docs/TASKGRAPH.md").write_text("\n".join(taskgraph_content), encoding="utf-8")
    print("✓ Generated docs/TASKGRAPH.md")

    # Generate DECISIONS.md
    decisions_content = generate_decisions()
    Path("docs/DECISIONS.md").write_text(decisions_content, encoding="utf-8")
    print("✓ Generated docs/DECISIONS.md")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
