use chematic::rxn::find_reaction_matches;
use xenosite_forest::mol::parse_mol;
use xenosite_forest::smarts::smarts_matches;

fn main() {
    let mol = parse_mol("Clc1ccccc1").unwrap();
    // SMARTS door (not reaction)
    for q in ["[#6]", "[#17]", "[#17]-[#6]", "[Cl]", "[c]", "[Cl]-[c]"] {
        let n = smarts_matches(&mol, q).map(|h| h.len()).unwrap_or(0);
        println!("SMARTS {q}: {n}");
    }
    // Reaction door
    for s in [
        "[#6:1]>>[#6:1]",
        "[c:1]>>[c:1]",
        "[C:1]>>[C:1]",
        "[#17:1]>>[#17:1]",
        "[Cl:1]>>[Cl:1]",
        "[#17:1]-[#6:2]>>[#17:1]-[#6:2]",
        "[Cl:1]-[c:2]>>[Cl:1]-[c:2]",
        "[Cl:1]-[C:2]>>[Cl:1]-[C:2]",
    ] {
        match find_reaction_matches(s, &[&mol]) {
            Ok(ms) => println!("RXN {s}: {}", ms.len()),
            Err(e) => println!("RXN {s}: ERR {e}"),
        }
    }
}
