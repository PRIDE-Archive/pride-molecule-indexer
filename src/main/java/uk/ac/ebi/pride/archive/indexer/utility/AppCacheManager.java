package uk.ac.ebi.pride.archive.indexer.utility;

import org.ehcache.Cache;
import org.ehcache.CacheManager;
import org.ehcache.config.Configuration;
import org.ehcache.config.builders.CacheManagerBuilder;
import org.ehcache.xml.XmlConfiguration;

import java.net.URL;
import java.util.List;

public class AppCacheManager {
    private static AppCacheManager instance = null;
    private CacheManager cacheManage = null;

    private AppCacheManager() {
        URL ehcacheXmlUrl = Thread.currentThread().getContextClassLoader().getResource("ehcache.xml");
        Configuration xmlConfig = new XmlConfiguration(ehcacheXmlUrl);
        cacheManage = CacheManagerBuilder.newCacheManager(xmlConfig);
        cacheManage.init();
        System.out.println("Cache Initialized");
    }

    public static AppCacheManager getInstance() {
        if (instance == null) {
            instance = new AppCacheManager();
        }
        return instance;
    }

    public Cache<Integer, Integer> getClustersCache(){
        return this.cacheManage.getCache("Clusters", Integer.class, Integer.class);
    }

    public Cache<String, Long> getPrideJsonSpectra() {
        return this.cacheManage.getCache("InputSpectra", String.class, Long.class);
    }


    public Cache<String, ? extends List> getProteinToPsmsCache() {
        return this.cacheManage.getCache("ProteinToPsms", String.class, List.class);
    }
}
