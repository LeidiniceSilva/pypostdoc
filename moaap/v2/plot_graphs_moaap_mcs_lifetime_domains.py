# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "March 03, 2026"
__description__ = "This script plot MCSs by lifetime duration"

import os
import pickle
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")


def open_mcs_era5(domain, start='2000-01', end='2009-12'):

    path = f'/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/ERA5/{domain}/output/'

    dates = pd.date_range(start=start, end=end, freq='MS')
    mcs = {}

    for d in dates:
        f = 'MCSs_' + d.strftime('%Y%m') + '__dt-1h_MOAAP-masks.pkl'
        f = os.path.join(path, f)

        if os.path.exists(f):
            with open(f, 'rb') as file:
                mcs[d.strftime('%Y-%m')] = pickle.load(file)

    return mcs


def open_mcs_cpm(domain, start='2000-01', end='2009-12'):

    path = f'/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/CPMs/{domain}/output/'

    dates = pd.date_range(start=start, end=end, freq='MS')
    mcs = {}

    for d in dates:
        f = 'MCSs_' + d.strftime('%Y%m') + '__dt-1h_MOAAP-masks.pkl'
        f = os.path.join(path, f)

        if os.path.exists(f):
            with open(f, 'rb') as file:
                mcs[d.strftime('%Y-%m')] = pickle.load(file)

    return mcs


def count_mcs_by_lifetime(mcs_charac):
    """
    Count number of MCSs by lifetime duration categories:
    5-10, 10-15, 15-20, >20 hours
    """
    # Initialize counters
    counts = {
        '5-10': 0,
        '10-15': 0,
        '15-20': 0,
        '>20': 0
    }
    
    for obj in mcs_charac.keys():
        for m in mcs_charac[obj].keys():
            # Get lifetime (duration) of each MCS
            # Assuming 'times' gives the time steps of the MCS
            times = mcs_charac[obj][m]['times']
            
            if len(times) > 0:
                # Calculate lifetime in hours
                lifetime_hours = len(times)
                
                # Categorize by lifetime
                if 5 <= lifetime_hours < 10:
                    counts['5-10'] += 1
                elif 10 <= lifetime_hours < 15:
                    counts['10-15'] += 1
                elif 15 <= lifetime_hours < 20:
                    counts['15-20'] += 1
                elif lifetime_hours >= 20:
                    counts['>20'] += 1
    
    return counts


def plot_duration_bars(ax, era5_counts, cpm_counts, domain_name, font_size=10):
    """
    Plot bar chart for MCS duration categories
    """
    categories = ['5-10', '10-15', '15-20', '>20']
    x = np.arange(len(categories))
    width = 0.35
    
    # Convert to list of values in the correct order
    era5_values = [era5_counts[cat] for cat in categories]
    cpm_values = [cpm_counts[cat] for cat in categories]
    
    # Create bars
    bars1 = ax.bar(x - width/2, era5_values, width, label='ERA5', 
                   color='red', alpha=0.75, edgecolor='red')
    bars2 = ax.bar(x + width/2, cpm_values, width, label='RegCM5', 
                   color='blue', alpha=0.75, edgecolor='blue')
    
    # Add value labels on top of bars
    for bar in bars1:
        height = bar.get_height()
        if height > 0:
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{int(height)}', ha='center', va='bottom', fontsize=8)
    
    for bar in bars2:
        height = bar.get_height()
        if height > 0:
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{int(height)}', ha='center', va='bottom', fontsize=8)
    
    # Customize plot
    ax.set_xlabel('Lifetime (hours)', fontsize=font_size, fontweight='bold')
    ax.set_ylabel('Number of MCSs', fontsize=font_size, fontweight='bold')
    ax.set_title(domain_name, loc='left', fontsize=font_size)
    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=font_size)
    ax.legend(loc=9, frameon=False, fontsize=font_size)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', alpha=0.3, linestyle='--')


# Main execution
if __name__ == "__main__":
    
    print("Loading MCS data...")
    
    # Load data for all domains
    domains = ['CAR-4', 'CSAM-3', 'EURR-3']
    
    era5_data = {}
    cpm_data = {}
    
    for domain in domains:
        print(f"  Loading {domain}...")
        era5_data[domain] = open_mcs_era5(domain)
        cpm_data[domain] = open_mcs_cpm(domain)
    
    # Count MCSs by lifetime for each domain
    print("\nCounting MCSs by lifetime duration...")
    
    era5_counts = {}
    cpm_counts = {}
    
    for domain in domains:
        era5_counts[domain] = count_mcs_by_lifetime(era5_data[domain])
        cpm_counts[domain] = count_mcs_by_lifetime(cpm_data[domain])
    
    # Print results
    print("\n" + "="*60)
    print("MCS COUNT BY LIFETIME DURATION (2000-2009)")
    print("="*60)
    
    for domain in domains:
        print(f"\n{domain}:")
        print(f"  {'Category':<10} {'ERA5':<10} {'RegCM5':<10}")
        print(f"  {'-'*10} {'-'*10} {'-'*10}")
        
        categories = ['5-10', '10-15', '15-20', '>20']
        for cat in categories:
            print(f"  {cat:<10} {era5_counts[domain][cat]:<10} {cpm_counts[domain][cat]:<10}")
        
        # Calculate totals
        era5_total = sum(era5_counts[domain].values())
        cpm_total = sum(cpm_counts[domain].values())
        print(f"  {'Total':<10} {era5_total:<10} {cpm_total:<10}")
    
    # Create figure with 1x3 subplots
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    font_size = 10
    
    # Plot for each domain
    for idx, domain in enumerate(domains):
        plot_duration_bars(axes[idx], era5_counts[domain], cpm_counts[domain], 
                          f'({chr(97+idx)}) {domain}', font_size)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure
    path_out = '/leonardo/home/userexternal/mdasilva/leonardo_work/MOAAP/figs'
    os.makedirs(path_out, exist_ok=True)
    name_out = f'pyplt_graphs_moaap_mcs_lifetime_distribution_2000-2009.png'
    plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
    
    print(f"\nFigure saved to: {os.path.join(path_out, name_out)}")
    plt.show()
    
    # Optional: Create a stacked bar chart for better comparison
    fig2, ax2 = plt.subplots(figsize=(10, 6))
    
    categories = ['5-10', '10-15', '15-20', '>20']
    x = np.arange(len(domains))
    width = 0.35
    
    # Prepare data for grouped bars by domain
    era5_by_domain = {cat: [era5_counts[d][cat] for d in domains] for cat in categories}
    cpm_by_domain = {cat: [cpm_counts[d][cat] for d in domains] for cat in categories}
    
    # Plot stacked bars
    bottom_era5 = np.zeros(len(domains))
    bottom_cpm = np.zeros(len(domains))
    
    colors = ['#ff9999', '#ffcccc', '#ff6666', '#ff3333']
    
    for i, cat in enumerate(categories):
        ax2.bar(x - width/2, era5_by_domain[cat], width, 
                label=f'ERA5 {cat}h', bottom=bottom_era5, color=colors[i], alpha=0.8)
        bottom_era5 += np.array(era5_by_domain[cat])
        
        ax2.bar(x + width/2, cpm_by_domain[cat], width,
                label=f'RegCM5 {cat}h', bottom=bottom_cpm, color=colors[i], alpha=0.5)
        bottom_cpm += np.array(cpm_by_domain[cat])
    
    ax2.set_xlabel('Domain', fontsize=font_size, fontweight='bold')
    ax2.set_ylabel('Number of MCSs', fontsize=font_size, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(domains, fontsize=font_size)
    ax2.legend(loc=9, ncol=2, fontsize=8)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.grid(axis='y', alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    
    name_out2 = f'pyplt_graphs_moaap_mcs_lifetime_domains_2000-2009.png'
    plt.savefig(os.path.join(path_out, name_out2), dpi=400, bbox_inches='tight')
    print(f"Stacked figure saved to: {os.path.join(path_out, name_out2)}")
    
    plt.show()
    exit()
