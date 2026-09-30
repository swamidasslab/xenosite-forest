use xenosite_forest::{
    accept_all_rules, accept_all_sites, as_forest_mol, atom_diff, phase_one,
};

fn main() {
    let set = phase_one();

    println!("=== N-OH-IQ ===");
    let r = as_forest_mol("Cn1c(=N)[nH]c2c3cccnc3ccc21").unwrap();
    let t = as_forest_mol("Cn1c(NO)nc2c3cccnc3ccc21").unwrap();
    let diff = atom_diff(r.mol(), t.mol());
    println!("R={} P={} root={}", r.csmi(), t.csmi(), diff.cost());
    for em in set.metabolize(&r, &accept_all_rules, &accept_all_sites, true) {
        let em = em.unwrap();
        if em.pattern_name != "hydroxylamine" {
            continue;
        }
        for p in &em.products {
            let d2 = atom_diff(p.mol(), t.mol());
            println!(
                "  hydroxylamine site={} → {} cost {}→{} hit={}",
                em.site,
                p.csmi(),
                diff.cost(),
                d2.cost(),
                p.csmi().as_ref() == t.csmi().as_ref()
            );
        }
    }

    println!("\n=== caffeine → theobromine (cleaving improvers) ===");
    let r = as_forest_mol("Cn1c(=O)c2c(ncn2C)n(C)c1=O").unwrap();
    let t = as_forest_mol("Cn1cnc2c(=O)[nH]c(=O)n(C)c12").unwrap();
    let diff = atom_diff(r.mol(), t.mol());
    println!("R={} P={} root={}", r.csmi(), t.csmi(), diff.cost());
    for em in set.metabolize(&r, &accept_all_rules, &accept_all_sites, true) {
        let em = em.unwrap();
        if !em.cleaves {
            continue;
        }
        for p in &em.products {
            let d2 = atom_diff(p.mol(), t.mol());
            if d2.cost() >= diff.cost() {
                continue;
            }
            println!(
                "  {} site={} → {} cost {}",
                em.pattern_name,
                em.site,
                p.csmi(),
                d2.cost()
            );
        }
    }
}
